from typing import Dict, Union

from desdeo_emo.EAs.BaseEA import BaseDecompositionEA, eaError
from desdeo_emo.population.Population import Population

#from desdeo_emo.selection.APD_Select import APD_Select #I uncommented this to test and commented the next line
from desdeo_emo.selection.APD_Select_constraints import APD_Select
from desdeo_emo.selection.Prob_APD_Select_MC import Prob_APD_select_MC  # superfast y considering MC samples
#from desdeo_emo.selection.Prob_Hybrid_APD_Select import Prob_Hybrid_APD_Select   # hybrid approach with mean selection
from desdeo_problem import MOProblem
import numpy as np
from desdeo_tools.interaction import (
    SimplePlotRequest,
    ReferencePointPreference,
    validate_ref_point_data_type,
    validate_ref_point_dimensions,
    validate_ref_point_with_ideal,
)

class RVEA(BaseDecompositionEA):
    """The python version reference vector guided evolutionary algorithm.

    Most of the relevant code is contained in the super class. This class just assigns
    the APD selection operator to BaseDecompositionEA.

    NOTE: The APD function had to be slightly modified to accomodate for the fact that
    this version of the algorithm is interactive, and does not have a set termination
    criteria. There is a time component in the APD penalty function formula of the type:
    (t/t_max)^alpha. As there is no set t_max, the formula has been changed. See below,
    the documentation for the argument: penalty_time_component

    See the details of RVEA in the following paper

    R. Cheng, Y. Jin, M. Olhofer and B. Sendhoff, A Reference Vector Guided
    Evolutionary Algorithm for Many-objective Optimization, IEEE Transactions on
    Evolutionary Computation, 2016

    Parameters
    ----------
    problem : MOProblem
        The problem class object specifying the details of the problem.
    population_size : int, optional
        The desired population size, by default None, which sets up a default value
        of population size depending upon the dimensionaly of the problem.
    population_params : Dict, optional
        The parameters for the population class, by default None. See
        desdeo_emo.population.Population for more details.
    initial_population : Population, optional
        An initial population class, by default None. Use this if you want to set up
        a specific starting population, such as when the output of one EA is to be
        used as the input of another.
    alpha : float, optional
        The alpha parameter in the APD selection mechanism. Read paper for details.
    lattice_resolution : int, optional
        The number of divisions along individual axes in the objective space to be
        used while creating the reference vector lattice by the simplex lattice
        design. By default None
    a_priori : bool, optional
        A bool variable defining whether a priori preference is to be used or not.
        By default False
    interact : bool, optional
        A bool variable defining whether interactive preference is to be used or
        not. By default False
    n_iterations : int, optional
        The total number of iterations to be run, by default 10. This is not a hard
        limit and is only used for an internal counter.
    n_gen_per_iter : int, optional
        The total number of generations in an iteration to be run, by default 100.
        This is not a hard limit and is only used for an internal counter.
    total_function_evaluations :int, optional
        Set an upper limit to the total number of function evaluations. When set to
        zero, this argument is ignored and other termination criteria are used.
    penalty_time_component: Union[str, float], optional
        The APD formula had to be slightly changed.
        If penalty_time_component is a float between [0, 1], (t/t_max) is replaced by
        that constant for the entire algorithm.
        If penalty_time_component is "original", the original intent of the paper is
        followed and (t/t_max) is calculated as
        (current generation count/total number of generations).
        If penalty_time_component is "function_count", (t/t_max) is calculated as
        (current function evaluation count/total number of function evaluations)
        If penalty_time_component is "interactive", (t/t_max)  is calculated as
        (Current gen count within an iteration/Total gen count within an iteration).
        Hence, time penalty is always zero at the beginning of each iteration, and one
        at the end of each iteration.
        Note: If the penalty_time_component ever exceeds one, the value one is used as
        the penalty_time_component.
        If no value is provided, an appropriate default is selected.
        If `interact` is true, penalty_time_component is "interactive" by default.
        If `interact` is false, but `total_function_evaluations` is provided,
        penalty_time_component is "function_count" by default.
        If `interact` is false, but `total_function_evaluations` is not provided,
        penalty_time_component is "original" by default.
    """

    def __init__(
        self,
        problem: MOProblem,
        population_size: int = None,
        population_params: Dict = None,
        initial_population: Population = None,
        alpha: float = 2,
        lattice_resolution: int = None,
        selection_type: str = None,
        interact: bool = False,
        use_surrogates: bool = False,
        n_iterations: int = 10,
        n_gen_per_iter: int = 100,
        total_function_evaluations: int = 0,
        time_penalty_component: Union[str, float] = None,
        keep_archive: bool = False,
        save_non_dominated: bool = False,
    ):
        super().__init__(
            problem=problem,
            population_size=population_size,
            population_params=population_params,
            initial_population=initial_population,
            lattice_resolution=lattice_resolution,
            interact=interact,
            use_surrogates=use_surrogates,
            n_iterations=n_iterations,
            n_gen_per_iter=n_gen_per_iter,
            total_function_evaluations=total_function_evaluations,
            keep_archive=keep_archive,
            save_non_dominated=save_non_dominated,
        )
        self.time_penalty_component = time_penalty_component
        time_penalty_component_options = ["original", "function_count", "interactive"]
        if time_penalty_component is None:
            if interact is True:
                time_penalty_component = "interactive"
            elif total_function_evaluations > 0:
                time_penalty_component = "function_count"
            else:
                time_penalty_component = "original"
        if not (type(time_penalty_component) is float or str):
            msg = (
                f"type(time_penalty_component) should be float or str"
                f"Provided type: {type(time_penalty_component)}"
            )
            eaError(msg)
        if type(time_penalty_component) is float:
            if (time_penalty_component <= 0) or (time_penalty_component >= 1):
                msg = (
                    f"time_penalty_component should either be a float in the range"
                    f"[0, 1], or one of {time_penalty_component_options}.\n"
                    f"Provided value = {time_penalty_component}"
                )
                eaError(msg)
            time_penalty_function = self._time_penalty_constant
        if type(time_penalty_component) is str:
            if time_penalty_component == "original":
                time_penalty_function = self._time_penalty_original
            elif time_penalty_component == "function_count":
                time_penalty_function = self._time_penalty_function_count
            elif time_penalty_component == "interactive":
                time_penalty_function = self._time_penalty_interactive
            else:
                msg = (
                    f"time_penalty_component should either be a float in the range"
                    f"[0, 1], or one of {time_penalty_component_options}.\n"
                    f"Provided value = {time_penalty_component}"
                )
                eaError(msg)
        self.time_penalty_function = time_penalty_function
        self.alpha = alpha
        self.selection_type = selection_type
        selection_operator = APD_Select(
            pop=self.population,
            time_penalty_function=self.time_penalty_function,
            alpha=alpha,
            selection_type=selection_type,
        )
        self.selection_operator = selection_operator

    def _time_penalty_constant(self):
        """Returns the constant time penalty value."""
        return self.time_penalty_component

    def _time_penalty_original(self):
        """Calculates the appropriate time penalty value, by the original formula."""
        return self._current_gen_count / self.total_gen_count

    def _time_penalty_interactive(self):
        """Calculates the appropriate time penalty value."""
        return self._gen_count_in_curr_iteration / self.n_gen_per_iter

    def _time_penalty_function_count(self):
        """Calculates the appropriate time penalty value."""
        return self._function_evaluation_count / self.total_function_evaluations



class ProbRVEA_MC(RVEA):
    def __init__(
        self,
        problem: MOProblem,
        population_size: int = None,
        population_params: Dict = None,
        initial_population: Population = None,
        alpha: float = 2,
        lattice_resolution: int = None,
        selection_type: str = None,
        interact: bool = False,
        use_surrogates: bool = False,
        n_iterations: int = 10,
        n_gen_per_iter: int = 100,
        total_function_evaluations: int = 0,
        time_penalty_component: Union[str, float] = None,
        keep_archive: bool = False,
        save_non_dominated: bool = False,
    ):
        super().__init__(
            problem=problem,
            population_size=population_size,
            population_params=population_params,
            initial_population=initial_population,
            lattice_resolution=lattice_resolution,
            interact=interact,
            use_surrogates=use_surrogates,
            n_iterations=n_iterations,
            n_gen_per_iter=n_gen_per_iter,
            total_function_evaluations=total_function_evaluations,
            keep_archive=keep_archive,
            save_non_dominated=save_non_dominated,
        )
        selection_operator = Prob_APD_select_MC(
            pop=self.population,
            time_penalty_function=self.time_penalty_function,
            alpha=alpha,
            selection_type=selection_type,
        )
        self.selection_operator = selection_operator

    """
    def _next_gen(self):
        size_pop = np.shape(self.population.objectives)[0]
        #print("Pop size before recombination: ",size_pop)
        if self.interact is True and size_pop <= self.thresh_size:
            self.population.individuals = self.iteration_archive_individuals
            self.population.objectives = self.iteration_archive_objectives
            self.population.fitness = self.iteration_archive_fitness
            self.population.uncertainity = self.iteration_archive_uncertainty
            print("population injected!!!!")
        else:
            offspring = self.population.mate()  # (params=self.params)
            self.population.add(offspring, self.use_surrogates)
            self._function_evaluation_count += offspring.shape[0]
        #print("Pop size after recombination: ",np.shape(self.population.objectives)[0])
        selected = self._select()
        #print("selected size:",np.shape(selected)[0])
        self.population.keep(selected)
        self._current_gen_count += 1
        self._gen_count_in_curr_iteration += 1
    """
"""
class HybridRVEA(RVEA):
    def __init__(
        self,
        problem: MOProblem,
        population_size: int = None,
        population_params: Dict = None,
        initial_population: Population = None,
        alpha: float = 2,
        lattice_resolution: int = None,
        a_priori: bool = False,
        interact: bool = False,
        use_surrogates: bool = False,
        n_iterations: int = 10,
        n_gen_per_iter: int = 100,
        total_function_evaluations: int = 0,
        time_penalty_component: Union[str, float] = None,
    ):
        super().__init__(
            problem=problem,
            population_size=population_size,
            population_params=population_params,
            initial_population=initial_population,
            lattice_resolution=lattice_resolution,
            a_priori=a_priori,
            interact=interact,
            use_surrogates=use_surrogates,
            n_iterations=n_iterations,
            n_gen_per_iter=n_gen_per_iter,
            total_function_evaluations=total_function_evaluations,
        )
        selection_operator = Prob_Hybrid_APD_Select(
            self.population, self.time_penalty_function, alpha
        )
        self.selection_operator = selection_operator
"""