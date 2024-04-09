from desdeo_problem.problem.Variable import Variable
from desdeo_problem.problem.Objective import ScalarObjective, VectorObjective
from desdeo_problem.problem.Problem import MOProblem
import math
import numpy as np
from Bio.SeqUtils import MeltingTemp

def primer_problem(s, US_GC_percent, US_Tm_f, US_Tm_r, Ct, US_hairpin, US_len_diff, US_Tm_diff) -> MOProblem:
    """An example on how to implement a problem with 3 objectives and 4 variables and no constraints.

    Returns:
        MOProblem: a problem object.
    """
    n_vars = 4
    deltaH = {
    "AA": -7.9,
    "TT": -7.9,
    "AT": -7.2,
    "TA": -7.2,
    "CA": -8.5,
    "GT": -8.4,
    "CT": -7.8,
    "GA": -8.2,
    "CG": -10.6,
    "GC": -9.8,
    "GG": -8.0,
    "CC": -8.0,
    "initGC": 0.1,
    "initAT": 2.3,
    "symmetryH": 0.0,
    }

    deltaS = {
        "AA": -22.2,
        "TT": -22.2,
        "AT": -20.4,
        "TA": -21.3,
        "CA": -22.7,
        "GT": -22.4,
        "CT": -21.0,
        "GA": -22.2,
        "CG": -27.2,
        "GC": -24.4,
        "GG": -19.9,
        "CC": -19.9,
        "initGC": -2.8,
        "initAT": 4.1,
        "symmetryS": -1.4,
    }

    bonds = {"A": "T", "T": "A", "G": "C", "C": "G"}
    rev_s = ""
    for i in s:
        rev_s = rev_s+ bonds[i]
    tl = len(s)
    pmin = 16
    pmax = 28
    plmin = 200
    plmax = 1000
    max_num= 10e10

    def Cal_delH_delS(s: str):
        delH = 0
        delS = 0
        nn = ""
        for i in range(1, len(s)):
            nn = s[i] + s[i + 1]
            delH += deltaH[nn]
            delS += deltaS[nn]
        delS += 0.368 * (len(s) - 1) * 50 * 0.001
        return delH, delS


    def TmSAN_old(s: str, Ct: float):
        dH, dS = Cal_delH_delS(s)
        return dH * 1000 / (dS + 1.987 * math.log(Ct / 4)) - 273.15
    
    def TmSAN(s,Ct):
        return MeltingTemp.Tm_NN(seq=s,dnac1=Ct,dnac2=0.0)
    
    def Tm_NN1(s,Ct):
        return MeltingTemp.Tm_NN(seq=s,dnac1=Ct,dnac2=0.0, nn_table=MeltingTemp.DNA_NN1)
    
    def Tm_NN2(s,Ct):
        return MeltingTemp.Tm_NN(seq=s,dnac1=Ct,dnac2=0.0, nn_table=MeltingTemp.DNA_NN2)
    
    def Tm_NN3(s,Ct):
        return MeltingTemp.Tm_NN(seq=s,dnac1=Ct,dnac2=0.0, nn_table=MeltingTemp.DNA_NN3)
    
    def Tm_NN4(s,Ct):
        return MeltingTemp.Tm_NN(seq=s,dnac1=Ct,dnac2=0.0, nn_table=MeltingTemp.DNA_NN4)
    
    def Tm_Wallace(s):
        return MeltingTemp.Tm_Wallace(seq=s)
    
    def Tm_GC(s):
        return MeltingTemp.Tm_GC(seq=s)

    def GC_percent(s: str):
        G_no = 0
        C_no = 0
        for c in s:
            if c == "G":
                G_no += 1
            elif c == "C":
                C_no += 1

        return ((G_no + C_no) / len(s))*100

    def selfdimer(s1):
        pairs = {('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')}
        s2 = s1[::-1]
        global extra
        extra = 0

        def solve(s1, s2, st):
            global extra
            extra = extra
            n = len(s1)
            m = len(s2)
            dimer = 0
            
            for i in range(st, n):
                cnt = 0
                res = 0
                
                for j in range(m):
                    if i + j >= n:
                        break
                    
                    if (s1[i + j], s2[j]) in pairs:
                        cnt += 1
                    else:
                        if cnt > 2:
                            res += 1
                        cnt = 0
                
                if cnt > 2:
                    res += 1
                
                if i == 0:
                    extra += res
                
                dimer += res
            
            return dimer

        a = solve(s1, s2[::-1], 0)
        b = solve(s2, s1[::-1], 1)

        c = extra // 2 + 1 if extra % 2 == 1 else extra // 2

        return (a + b - c)
    
    def crossdimer(s1, s2):
        pairs = {('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')}
        global extra
        extra = 0

        def solve(s1, s2, st):
            global extra
            extra = extra
            n = len(s1)
            m = len(s2)
            dimer = 0
            
            for i in range(st, n):
                cnt = 0
                res = 0
                
                for j in range(m):
                    if i + j >= n:
                        break
                    
                    if (s1[i + j], s2[j]) in pairs:
                        cnt += 1
                    else:
                        if cnt > 2:
                            res += 1
                        cnt = 0
                
                if cnt > 2:
                    res += 1
                
                if i == 0:
                    extra += res

                # print(res)
                dimer += res
            
            return dimer

        a = solve(s1, s2[::-1], 0)
        b = solve(s2, s1[::-1], 0)

        return (a + b)

    def pair(s1, s2):
        pairs = {('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')}
        
        n = len(s1)
        m = len(s2)
        dimer = 0
        
        for i in range(1, n):
            cnt = 0
            res = 0
            
            for j in range(m):
                if i + j >= n:
                    break
                
                if (s1[i + j], s2[j]) in pairs:
                    cnt += 1
                else:
                    if cnt == m:
                        res += 1
                    cnt = 0
            
            if cnt == m:
                res += 1
            
            dimer += res
        
        return dimer > 0

    def hairpin(s):
        n = len(s)
        ans = 0
        length = 0
        
        for j in range(n // 2 + 1):
            val = 0
            
            for i in range(2, n):
                s1 = s[j: j + max(length, i + 1)]
                s2 = ""
                
                if j + i + 1 < n:
                    s2 = s[j + i + 1:]
                
                s1 = s1[::-1]
                
                if pair(s2, s1):
                    val += 1
                    length = max(length, i + 1)
            
            if val > 0:
                ans += 1
        if ans <= US_hairpin:
            ans = 0
        else:
            ans = ans - US_hairpin
        return ans



    def fprimer(x):
        #rs = x[0] + x[2] - x[3]
        #return s[rs : (rs + x[3])]
        #print(int(x[0]))
        #print(int(x[0]+x[1]))
        #print(s[int(x[0]):int(x[0]+x[1])])
        #return rev_s[int(x[0]):int(x[0]+x[1])]
        return s[int(x[0]):int(x[0]+x[1])]

    def rprimer(x):
        #return rev_s[x[0] : (x[0] + x[1] + 1)]
        #print(int(x[0]+x[2]-x[3]))
        #print(int(x[0]+x[2]))
        #print(s[int(x[0]+x[2]-x[3]): int(x[0]+x[3])])
        #return s[int(x[0]+x[2]-x[3]): int(x[0]+x[2])]
        return rev_s[int(x[0]+x[2]-x[3]): int(x[0]+x[2])][::-1]

    def Tm_all(s, Ct):
        return np.asarray([Tm_NN1(s, Ct),Tm_NN2(s,Ct),Tm_NN3(s,Ct), Tm_NN4(s,Ct), Tm_GC(s), Tm_Wallace(s)])

    def f_1(fprimer):
        #return abs(TmSAN(fprimer, Ct) - US_Tm), 0
        tm_all = np.abs(Tm_all(fprimer,Ct)-US_Tm_f)
        return np.mean(tm_all), np.std(tm_all)


    def f_2(rprimer):
        #return abs(TmSAN(rprimer, Ct) - US_Tm), 0
        tm_all = np.abs(Tm_all(rprimer,Ct)-US_Tm_r)
        return np.mean(tm_all), np.std(tm_all)

    def f_3(fprimer, rprimer):
        """
        diff_Tm = abs(TmSAN(fprimer, Ct) - TmSAN(rprimer, Ct))
        if diff_Tm <= US_Tm_diff:
            return 0, 0
        return diff_Tm - US_Tm_diff, 0
        """
        tm_all_f = Tm_all(fprimer,Ct)
        tm_all_r = Tm_all(rprimer,Ct)
        diff_Tm = np.abs(tm_all_f-tm_all_r) - US_Tm_diff
        return np.mean(diff_Tm), np.std(diff_Tm)



    def f_4(fprimer, rprimer):
        # GCclamp

        if fprimer[-1] in "GC" and rprimer[-1] in "GC":
            GCclamp = 0
        elif fprimer[-1] in "GC" or rprimer[-1] in "GC":
            GCclamp = 1
        else:
            GCclamp = 2

        # len_diff
        diff_len = abs(len(fprimer) - len(rprimer))
        if diff_len <= US_len_diff:
            len_diff = 0
        else:
            len_diff = diff_len - US_len_diff

        return GCclamp + len_diff, 0

    def f_5(fprimer):
        return abs(GC_percent(fprimer) - US_GC_percent), 0


    def f_6(rprimer):
        return abs(GC_percent(rprimer) - US_GC_percent), 0

    def f_7(fprimer, rprimer):
        return selfdimer(fprimer) + selfdimer(rprimer) + crossdimer(fprimer, rprimer), 0

    def f_8(fprimer, rprimer):
        return hairpin(fprimer) + hairpin(rprimer), 0

    def constraint_check(x):
        if int(x[0]+x[1])>tl-1 or int(x[0]+x[2])>tl-1 or int(x[0])>int(x[0]+x[2]-x[3]):
            return False
        else:
            return True


    def f1(x: np.ndarray) -> np.ndarray:
        if constraint_check(x):
            return f_1(fprimer(x))
        else:
            return max_num, 0
        

    def f2(x: np.ndarray) -> np.ndarray:
        if constraint_check(x):
            return f_2(rprimer(x))
        else:
            return max_num, 0
        

    def f3(x: np.ndarray) -> np.ndarray:
        if constraint_check(x):
            return f_3(fprimer(x), rprimer(x))
        else:
            return max_num, 0
    
    def f4(x: np.ndarray) -> np.ndarray:
        if constraint_check(x):
            return f_4(fprimer(x), rprimer(x))
        else:
            return max_num, 0
    
    def f5(x: np.ndarray) -> np.ndarray:
        if constraint_check(x):
            return f_5(fprimer(x))
        else:
            return max_num, 0
           
    def f6(x: np.ndarray) -> np.ndarray:
        if constraint_check(x):
            return f_6(rprimer(x))
        else:
            return max_num, 0

    def f7(x: np.ndarray) -> np.ndarray:
        if constraint_check(x):
            return f_7(fprimer(x), rprimer(x))
        else:
            return max_num, 0

    def f8(x: np.ndarray) -> np.ndarray:
        if constraint_check(x):
            return f_8(fprimer(x), rprimer(x))
        else:
            return max_num, 0
        
    def objective_modified(x):
        #if constraint_check(x):
        return f1(x),f2(x),f3(x),f4(x),f5(x),f6(x),f7(x),f8(x)
        #else:
        #    return np.inf()
        
    def vect_f(x):
    #x=np.random.logistic(1-0.5*(1-x))
        if isinstance(x, list):
            if len(x) == n_vars:
                return [objective_modified(x)]
            elif len(x[0]) == n_vars:
                return list(map(objective_modified, x))
        else:
            if x.ndim == 1:
                return [objective_modified(x)]
            elif x.ndim == 2:
                return list(map(objective_modified, x))
        raise TypeError("Unforseen problem, contact developer")
    
    
    """
    # args: name of objective, evaluator
    objective_1 = ScalarObjective("F1", f1)
    objective_2 = ScalarObjective("F2", f2)
    objective_3 = ScalarObjective("F3", f3)
    objective_4 = ScalarObjective("F4", f4)
    objective_5 = ScalarObjective("F5", f5)
    objective_6 = ScalarObjective("F6", f6)
    
    objectives = [objective_1, objective_2, objective_3, objective_4, objective_5, objective_6]
    """
    #vec_obj = VectorObjective(["f1,f2,f3,f4,f5,f6"])
    # args: name of variable, initial value, lower bound, upper bound
    variable_1 = Variable("x1", 0, 0, tl-pmin)
    variable_2 = Variable("x2", pmin, pmin, pmax)
    variable_3 = Variable("x3", 500, 500, 800)
    variable_4 = Variable("x4", pmin, pmin, pmax)

    variables = [variable_1, variable_2, variable_3, variable_4]
    obj_names = ["f" + str(i + 1) for i in range(8)]
    objectives_vect = VectorObjective(name=obj_names, evaluator=vect_f)

    # instantiate MOProblem object
    problem = MOProblem(variables=variables, objectives=[objectives_vect])
    #problem = MOProblem(vec_obj, variables)

    return problem

    