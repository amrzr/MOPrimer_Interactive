import numpy as np
from warnings import warn
from typing import List, Callable
from desdeo_emo.selection.SelectionBase import InteractiveDecompositionSelectionBase
from desdeo_emo.population.Population import Population
from typing import TYPE_CHECKING
from desdeo_emo.utilities.ProbabilityWrong import Probability_wrong
import os
os.environ["OMP_NUM_THREADS"] = "1"


class Prob_APD_select_MC(InteractiveDecompositionSelectionBase):
    """The selection operator for the RVEA algorithm. Read the following paper for more
        details.
        R. Cheng, Y. Jin, M. Olhofer and B. Sendhoff, A Reference Vector Guided
        Evolutionary Algorithm for Many-objective Optimization, IEEE Transactions on
        Evolutionary Computation, 2016
    Parameters
    ----------
    pop : Population
        The population instance
    time_penalty_function : Callable
        A function that returns the time component in the penalty function.
    alpha : float, optional
        The RVEA alpha parameter, by default 2
    """

    def __init__(
        self, pop: Population,
        time_penalty_function: Callable,
        alpha: float = 2,
        selection_type: str = None,
    ):
        super().__init__(pop.pop_size, pop.problem.n_of_fitnesses, selection_type)
        self.time_penalty_function = time_penalty_function
        self.alpha = alpha
        self.n_of_objectives = pop.problem.n_of_objectives

    def do(self, pop: Population) -> List[int]:
        
        self.request_preferences
        fitness = pop.fitness
        uncertainty = pop.uncertainity
        penalty_factor = self._partial_penalty_factor()
        refV = self.vectors.neighbouring_angles_current
        fmin = np.amin(fitness, axis=0)
        translated_fitness = fitness - fmin
        if self.pref_refpnt is None:
            translated_refpnt = np.zeros_like(fmin)
        else:
            translated_refpnt = self.pref_refpnt - fmin
        pwrong = Probability_wrong(mean_values=translated_fitness, stddev_values=uncertainty, n_samples=1000)
        pwrong.vect_sample_f()

        fitness_norm = np.linalg.norm(pwrong.f_samples, axis=1)
        fitness_norm = np.repeat(np.reshape(fitness_norm, (len(fitness), 1, pwrong.n_samples)), len(fitness[0, :]), axis=1)

        normalized_fitness = np.divide(pwrong.f_samples, fitness_norm)  # Checked, works.


        #Find cosine angles for all the samples
        cosine = np.tensordot(normalized_fitness, np.transpose(self.vectors.values), axes=([1], [0]))
        cosine = np.transpose(cosine,(0,2,1))

        if cosine[np.where(cosine > 1)].size:
            cosine[np.where(cosine > 1)] = 1
        if cosine[np.where(cosine < 0)].size:
            cosine[np.where(cosine < 0)] = 0
        # Calculation of angles between reference vectors and solutions
        theta = np.arccos(cosine)
        # Reference vector asub_population_indexssignment
        #pwrong.compute_pdf(cosine)
        # Compute rank of cos theta (to be vectorized)
        rank_cosine = np.mean(cosine,axis=2)
        assigned_vectors = np.argmax(rank_cosine, axis=1)
        selection = np.array([], dtype=int)

        vector_selection = None

        for i in range(0, len(self.vectors.values)):
            sub_population_index = np.atleast_1d(
                np.squeeze(np.where(assigned_vectors == i))
            )
            sub_population_fitness = pwrong.f_samples[sub_population_index]

            if len(sub_population_fitness > 0):
                # APD Calculation
                angles = theta[sub_population_index, i]
                angles = np.divide(angles, refV[i])  # This is correct.
                # You have done this calculation before. Check with fitness_norm
                # Remove this horrible line
                
                #sub_pop_fitness_magnitude = np.sqrt(
                #    np.sum(np.power(sub_population_fitness, 2), axis=1)
                #)
                zz=np.repeat(translated_refpnt.reshape(1,-1),sub_population_fitness.shape[0], axis=0)
                zz2=np.repeat(np.array([zz]),sub_population_fitness.shape[2], axis=0).transpose((1,2,0))
                sub_pop_fitness_magnitude = np.sqrt(
                    np.sum(np.power((sub_population_fitness-zz2), 2), axis=1)
                )
                sub_popfm = np.reshape(sub_pop_fitness_magnitude, (1, len(sub_pop_fitness_magnitude[:,0]), pwrong.n_samples))
                angles = np.reshape(angles,(1,len(angles),pwrong.n_samples))


                #### Overall Mean/Median of apd or using MC samples
                apd = np.multiply(
                    sub_popfm,
                    (1 + np.dot(penalty_factor, angles))
                )
                #rank_apd = np.mean(apd, axis=2)
                rank_apd = pwrong.compute_rank_MC(apd)
                minidx = np.where(rank_apd[0] == np.nanmin(rank_apd[0]))

                if np.isnan(apd).all():
                    continue
                selx = sub_population_index[minidx]
                if selection.shape[0] == 0:
                    selection = np.hstack((selection, np.transpose(selx[0])))
                    vector_selection = np.asarray(i)
                else:
                    selection = np.vstack((selection, np.transpose(selx[0])))
                    vector_selection = np.hstack((vector_selection, i))

        if selection.shape[0] == 1:
            print("Only one individual!!")
            rand_select = np.random.randint(len(fitness), size=1)
            selection = np.vstack((selection, np.transpose(rand_select[0])))
        return selection.squeeze()


    def _partial_penalty_factor(self) -> float:
        """Calculate and return the partial penalty factor for APD calculation.
            This calculation does not include the angle related terms, hence the name.
            If the calculated penalty is outside [0, 1], it will round it up/down to 0/1

        Returns
        -------
        float
            The partial penalty value
        """
        
        if self.time_penalty_function() < 0:
            px = 0
        elif self.time_penalty_function() > 1:
            px = 1
        else:
            px= self.time_penalty_function()
        penalty = ((px) ** self.alpha) * self.n_of_objectives

        return penalty
        """
        penalty = ((self.time_penalty_function()) ** self.alpha) * self.n_of_objectives
        if penalty < 0:
            penalty = 0
        if penalty > 1:
            penalty = 1
        return penalty
        """

