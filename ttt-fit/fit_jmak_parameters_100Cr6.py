"""
Fit of temperature-dependent JMAK parameters N(T) and b(T) based on TTT diagram.

Three different variants of fitting are compared:

- The isothermal fit is performed for each individual temperature. 
  Because there are only 2 data point per temperature, and we fit 2 constant
  parameters, we should obtain a perfect match of the reference data.
  
- The parameter fit approximates the temperature-dependency of the 
  parameters by a polynomial ansatz for logb and N. The polynomial coefficients
  are optimised to reproduce the parameters obtained from the isothermal
  fit as best as possible.
  
- The model fit ignores the isothermal fit, and instead optimises the
  polynomial coefficients to minimize the differences in phase fractions
  after a given transformation time at each temperature.
  

"""
import numpy as np
from numpy.polynomial import Polynomial
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt

# comment this line if LaTeX is not installed
plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = "\n".join([r"\usepackage[utf8]{inputenc}",
                                                 r"\usepackage[T1]{fontenc}",
                                                 r'\usepackage{siunitx}'])


def read_csv_file(file_path, beta_vals):
    """
    Read a csv file and return a data frame .
    
    The csv file must contain two columns named 'Temperature' and 
    'log10(time)', with semicolon as a separator and comma as the decimal
    separator for float values.

    Parameters
    ----------
    file_path : str
        Path to the csv file.
    beta_vals : Sequence[float]
        Bainite percentage for each row with the same temperature, in 
        decreasing order.

    Returns
    -------
    df : DataFrame

    """
    # read csv data from image digitizer
    df = pd.read_csv(file_path,
                     decimal=',',
                     delimiter=';',
                     header=0)
    
    # make sure identical temperature values can be recognized
    df = df.round({"Temperature": 0})
    
    temps = df["Temperature"].unique()
    
    # convert times from log scale
    df["Time"] = 10**df["log10(time)"]
    
    # fit is performed in the unit system of the simulations using ms
    df["Time"] = 10**3*df["Time"]
    df = df.drop("log10(time)", axis=1)
    
    
    # add a new column for the percentage of bainite
    df["beta"] = np.nan
    
    for temp in temps:
        # get all indices for this temperature
        indices = np.flatnonzero(df["Temperature"] == temp)
        
        # bainite percentages for this temperature
        # assumes that highest values are missing (if any)
        percentages = beta_vals[-len(indices):]
        
        df.iloc[indices,2] = percentages
        
    return df


def jmak_beta(t, b, N):
    """
    Isothermal JMAK model solved for bainite fraction after given time t.
    
    Assumes that the transformation starts at t=0.
    
    Parameters
    ----------
    t : float
		Time since transformation started.
	b, N : float
		JMAK model parameters.
		
	Returns
	-------
	float
		Bainite volume fraction.
		
    """
    return 1 - np.exp(-b*t**N)


def jmak_t(beta, b, N):
    """
    Isothermal JMAK model solved for the time to reach given bainite fraction.
    
    Parameters
    ----------
    beta : float
		Volume fraction of bainite.
	b, N : float
		JMAK model parameters.
		
	Returns
	-------
	float
		Time since transformation started.
    
    """
    return (1/b*np.log(1/(1-beta)))**(1/N)

def logb_N_param(T, *params):
    """
    Polynomial ansatz for parameters logb(T) and N(T).
    
    Coefficients can be for any order of polynomial and must be sorted 
    (logb0, logb1, .. ,N0, N1, ...).
    
    Parameters
    ----------
    T : float
		Temperature.
    params : Sequence[float]
		Polynomial coefficients.
	
	Returns
	-------
	logb, N : float
		JMAK model parameters logb and N evaluated at temperature T.
    
    """
    # infer the order of the polynomial
    order = int(len(params)/2) # 1 for constant, etc
    
    # extract coefficients for logb and N
    params_logb = np.array(params[0:order])
    params_N = np.array(params[order:])
    
    # powers of temperature
    Tpot = np.array([T**i for i in range(order)])
    
    # evaluate parametrisation
    logb = params_logb @ Tpot
    N = params_N @ Tpot
    
    return logb, N

def jmak_beta_param(t, T, *params):
    """
    Isothermal, temperature-dependent JMAK model solved for bainite fraction.
    
    Model parameters logb and N are obtained from polynomial ansatz.
    
    Coefficients can be for any order of polynomial and must be sorted 
    (logb0, logb1, .. ,N0, N1, ...).
    
    Parameters
    ----------
    t : float
		Time since transformation started.
	T : float
		Temperature.
	params : Sequence[float]
		Polynomial coefficients.
		
	Returns
	-------
	float
		Bainite volume fraction.
		
    """
    logb, N = logb_N_param(T, *params)
    
    return jmak_beta(t, np.exp(logb), N)

def jmak_t_param(beta, T, *params):
    """
    Isothermal, temperature-dependent JMAK model solved for the time.
    
    Model parameters logb and N are obtained from polynomial ansatz.
    
    Coefficients can be for any order of polynomial and must be sorted 
    (logb0, logb1, .. ,N0, N1, ...).
    
    Parameters
    ----------
    beta : float
		Volume fraction of bainite.
	T : float
		Temperature.
	params : Sequence[float]
		Polynomial coefficients.
		
	Returns
	-------
	float
		Time since transformation started.
		
    """
    logb, N = logb_N_param(T, *params)
    
    return jmak_t(beta, np.exp(logb), N)


def logb_N_isothermal_fit(df):
    """
    Calculate the isothermal fit for discrete JMAK parameters logb_i and N_i.
    
    The isothermal fit is performed for each unique temperature in the
    DataFrame df individually and should recover the original data up to 
    numerical tolerances.
    
    Parameters
    ----------
    df : DataFrame
		Reference data. The DataFrame must have columns 'Temperature',
		'Time', and 'beta', with two data points for each temperature
		representing 1% and 99% transformation.
		
	Returns
	-------
	T : ndarray
		Unique temperatures T_i.
	logb : ndarray
		JMAK model parameter logb_i at each temperature T_i.
	N : ndarray
		JMAK model parameter N_i at each temperature T_i.
		    
    """
    # the data at different temperatures is best reproduced with different
    # initial guesses, therefore we provide multiple initial guesses which
    # will be used at each temperature. The best fit is then extracted
    # from the results.
    logb0vec = (-6, -3.55)
    N0vec = (2.5, 1.74)
        
    # define bounds for parameters log(B) and N -  use +/- np.inf
    # for no bound. These values are arbitrarily chosen and increase the
    # numerical stability. If bounds are used, it should be checked that
    # the model parameters actually reproduce the reference data.
    logb_min = -np.inf
    N_min = 1.5
    logb_max = np.inf
    N_max = 3.5
    
    bounds = ((logb_min, N_min),(logb_max, N_max))
    
    # get an array of unique temperature values in the data
    unique_temps = df["Temperature"].unique()
    
    # we model the function ln(1-beta(t)), dependent on parameters logb and N.
    # the model function for scipy.optimize.curve_fit must have the signature
    # f(t, *params)

    # models ln(1-beta)
    def model(t, logb, N):
        return -np.exp(logb)*t**N

    # initialize lists for optimal parameters at each temperature
    logb_list = []
    N_list = []

    # loop over unique temperatures and fit logb, N for each
    for temp in unique_temps:
        # get the data points at the given tempretaure - should be 2 data
        # points for each temp
        temp_df = df[df["Temperature"]==temp]
        time = temp_df["Time"]
        beta = temp_df["beta"]
        
        # loop over initial guesses and store optimization results
        results = []
        for (logb0, N0) in zip(logb0vec, N0vec):
            # optimize the parameters
            (logb, N), _ = curve_fit(model,
                                     time,
                                     np.log(1-beta),
                                     p0=(logb0, N0),
                                     bounds=bounds)
            
            results.append((logb, N))
            
        # cost for each result
        result_errors = [np.linalg.norm(beta - jmak_beta(time, np.exp(logb), N))
                         for logb, N in results]
        
        # extract best result and append it to optimal parameter lists
        imin = np.argmin(result_errors)
        logb, N = results[imin]
        
        logb_list.append(logb)
        N_list.append(N)
        
    return unique_temps, np.array(logb_list), np.array(N_list)


def logb_N_parameter_fit(T_i, logb_i, N_i, order):
    """
    Calculate the parameter fit for parameters p of logb(T) and N(T).
    
    The parameter fit minimises the deviation of the continuous polynomials
    logb_p(T) and N_p(T) from the discrete parameters logb_i and N_i.
    The discrete parameters should have been obtained prior to the
    function call by optimising at each unique temperature.
    
    Parameters
    ----------
    T_i : Sequence[float]
		Temperatures at which the reference parameters were obtained.
	logb_i, N_i : Sequence[float]
		Reference parameters at temperatures T_i.
	order : int
		Order of the polynomials used to approximate logb(T) and N(T).
		
	Returns
	-------
	logb_p, N_p : np.Polynomial
		Polynomial objects representing the parameter fit.
    
    """
    # the convert() method expands the polynomial form to match the
    # expressions given in the paper (to make coefficients match)
    logb_p = Polynomial.fit(T_i, logb_i, order).convert()
    N_p = Polynomial.fit(T_i, N_i, order).convert()

    return logb_p, N_p


def logb_N_model_fit(df, p0):
    """
    Calculate the model fit for parameters p of logb(T) and N(T).
    
    The model fit minimizes the deviation from 1% and 99% bainite at the
    corresponding times in the ITT diagram. It takes into account all data
    points in the DataFrame df at once and relies on a good initial guess p0.
        
    Parameters
    ----------
    df : DataFrame
		Reference data. The DataFrame must have columns 'Temperature',
		'Time', and 'beta', with two data points for each temperature
		representing 1% and 99% transformation.
    p0 : Sequence[float]
		Initial guess for polynomial coefficients of the form
		(p_logb0, p_logb1, .. ,p_N0, p_N1, ...). Polynomial order is
		inferred from the length of the initial guess.
    
    Returns
    -------
    logb_p, N_p : np.Polynomial
		Polynomial objects representing the model fit.
    
    """
    # rearrange the data so that x consists of tuples of time and
    # temperature, while y is the corresponding bainite fraction
    temperature = df["Temperature"]
    time = df["Time"]
    beta = df["beta"]
    
    xdata = np.vstack((time, temperature))
    
    # model function 
    def model(x, *params):
        t, T = x
        return jmak_beta_param(t, T, *params)
    
    # actual optimisation
    params, cov = curve_fit(model,
                            xdata,
                            beta,
                            p0=p0,
                            )
    
    # extract the polynomial coefficients for logb(T) and N(T) from the
    # optimisation result
    params_logb, params_N = params.reshape(2,-1)
    
    # create Polynomial objects
    logb_p = Polynomial(params_logb).convert()
    N_p = Polynomial(params_N).convert()
    
    return logb_p, N_p


def plot_parameters(T_i, logb_i, N_i,
                    T_param, logb_param, N_param,
                    T_model, logb_model, N_model,
                    ):
    """
    Create a figure with 2 subplots for parameters over temperature.
    
    The figure is saved in the current directory as 
    `plot_jmak_parameters.png`.
    
    Parameters
    ----------
    T_i, logb_i, N_i : Sequence[float]
		Discrete model parameters log_i and N_i and temperatures T_i.
	T_param, logb_param, N_param : Sequence[float]
		Discrete model parameters logb_param and N_param obtained from
		parameter fit at temperatures T_param. Since a polynomial
		ansatz is used, these are only discrete for the sake of plotting.
    T_model, logb_model, N_model : Sequence[float]
		Discrete model parameters logb_param and N_param obtained from
		model fit at temperatures T_param. Since a polynomial
		ansatz is used, these are only discrete for the sake of plotting.
		
    Returns
    -------
    None
    
    """
    _, (ax1, ax2) = plt.subplots(1, 2, figsize=(14/2.45,8/2.45))
    
    ax1.plot(T_i, logb_i, 'x', color="k", label=r"isothermal fit $\hat{b}_i$")
    ax1.plot(T_param, logb_param, color="tab:blue", label=r"parameter fit $\tilde{b}_p$")
    ax1.plot(T_model, logb_model, '--', color="tab:orange", label=r"model fit $b_p^*$")
    
    plt.sca(ax1)

    plt.xlabel(r"Temperature $T$ / \si{\celsius}")
    plt.ylabel(r"$\ln(b)$")
    plt.grid()
    plt.legend()
    
    ax2.plot(T_i, N_i, 'x', color="k", label=r"isothermal fit $\hat{N}_i$")
    ax2.plot(T_param, N_param, color="tab:blue", label=r"parameter fit $\tilde{N}_p$")
    ax2.plot(T_model, N_model, '--', color="tab:orange", label=r"model fit $N_p^*$")
    
    plt.sca(ax2)
    plt.xlabel(r"Temperature $T$ / \si{\celsius}")
    plt.ylabel("$N$")
    plt.grid()
    plt.legend()
    
    plt.tight_layout()
    
    plt.savefig("plot_jmak_parameters.png", dpi=600)
    
    
def plot_ttt_diagram_all_models(T_i,
                                t_ref_1, 
                                t_ref_99,
                                t_i_1,
                                t_i_99,
                                T_param,
                                t_param_1,
                                t_param_99,
                                T_model,
                                t_model_1,
                                t_model_99,
                                ):
    """
    Create a TTT diagram with reference data and both fits.
    
    The figure is saved in the current directory as 
    `plot_ttt_diagram_both_fits.png`.
    
    Parameters
    ----------
    T_i : Sequence[float]
		Discrete temperature values in the reference data.
	t_ref_1, t_ref_99 : Sequence[float]
		Transformation time for 1% and 99% bainite in the reference data.
	t_i_1, t_i_99 : Sequence[float]
		Transformation time for 1% and 99% bainite obtained by fitting
		the model parameters logb and N at each temperature in T_i
		separately (Verification of discrete reference parameters).
	T_param, t_param_1, t_param_99 : 
		Temperatures and transformation times for 1% and 99% bainite from
		the parameter fit. Since a polynomial ansatz is used, these are
		only discrete for the sake of plotting.
	T_model, t_model_1, t_model_99 : 
		Temperatures and transformation times for 1% and 99% bainite from
		the model fit. Since a polynomial ansatz is used, these are
		only discrete for the sake of plotting.
		
	Returns
	-------
	None
    
    """
    plt.figure(figsize=(10/2.45,8/2.45))
    
    # plot reference
    plt.plot(t_ref_1,
             T_i,
             '.',
             label="Reference",
             color="tab:gray",
             linewidth=1,
             )

    plt.plot(t_ref_99,
             T_i,
             '.',
             color="tab:gray",
              linewidth=2,
             )
    
    plt.plot(t_i_1,
             T_i,
             'x',
             label="isothermal fit",
             color="k")
    
    plt.plot(t_i_99,
             T_i,
             'x',
             color="k")
    
    plt.plot(t_param_1,
             T_param,
             label="parameter fit",
             color="tab:blue")
    
    plt.plot(t_param_99,
             T_param,
             color="tab:blue")
    
    plt.plot(t_model_1,
             T_model,
             '--',
             label="model fit",
             color="tab:orange")
    
    plt.plot(t_model_99,
             T_model,
             '--',
             color="tab:orange")
    
    
    plt.xlim([1, plt.xlim()[1]])
    plt.xscale("log")
    plt.xlabel(r"Time $t$ / \si{\second}")
    plt.ylabel(r"Temperature $T$  / \si{\celsius}")
    plt.legend()
    plt.grid()
    
    plt.tight_layout()
    
    plt.savefig("plot_ttt_diagram_both_fits.png", dpi=600)
    

def plot_ttt_diagram_model_fit(T_i,
                               t_ref_1, 
                               t_ref_99,
                               t_i_1,
                               t_i_99,
                               T_model,
                               t_model_1,
                               t_model_99,
                               ):
    """
    Create a TTT diagram with reference data and the model fit.
    
    The figure is saved in the current directory as 
    `plot_ttt_diagram_model_fit.png`.
    
    Parameters
    ----------
    T_i : Sequence[float]
		Discrete temperature values in the reference data.
	t_ref_1, t_ref_99 : Sequence[float]
		Transformation time for 1% and 99% bainite in the reference data.
	t_i_1, t_i_99 : Sequence[float]
		Transformation time for 1% and 99% bainite obtained by fitting
		the model parameters logb and N at each temperature in T_i
		separately (Verification of discrete reference parameters).
	T_model, t_model_1, t_model_99 : 
		Temperatures and transformation times for 1% and 99% bainite from
		the model fit. Since a polynomial ansatz is used, these are
		only discrete for the sake of plotting.
		
	Returns
	-------
	None
    
    
    """
    plt.figure(figsize=(10/2.45,8/2.45))
    
    # plot reference
    plt.plot(t_ref_1,
             T_i,
             '.',
             label=r"Reference \SI{1}{\percent} bainite",
             color="tab:gray",
             linewidth=1,
             )

    plt.plot(t_ref_99,
             T_i,
             'x',
             color="tab:gray",
             label=r"Reference \SI{99}{\percent} bainite",
             linewidth=2,
             )
    
    plt.plot(t_model_1,
             T_model,
             '--',
             label=r"Model \SI{1}{\percent} bainite",
             color="tab:orange")
    
    plt.plot(t_model_99,
             T_model,
             '-.',
             label=r"Model \SI{99}{\percent} bainite",
             color="tab:orange")
    
    
    plt.xlim([1, plt.xlim()[1]])
    plt.xscale("log")
    plt.xlabel(r"Time $t$ / \si{\second}")
    plt.ylabel(r"Temperature $T$  / \si{\celsius}")
    plt.legend()
    plt.grid()
    
    plt.tight_layout()
    
    plt.savefig("plot_ttt_diagram_model_fit.png", dpi=600)
    
    
def print_coefficients(logb_param, N_param, logb_model, N_model):
    """
    Print the coefficients of the optimized polynomials.
	
	Parameters
	----------
	logb_param, N_param : np.Polynomial
		Approximation of the model parameters from parameter fit.
	logb_model, N_model : np.polynomial
		Approximation of the model parameters from model fit.

	Returns
	-------
	None
	
    """
    print(f"logb_param: {logb_param}\n")
    print(f"N_param: {N_param}\n")
    print(f"logb_model: {logb_model}\n")
    print(f"N_model: {N_model}\n")
    
    logb_param_str = [f"{c:0.4e}" for c in logb_param.convert().coef]
    N_param_str = [f"{c:0.4e}" for c in N_param.convert().coef]
    logb_model_str = [f"{c:0.4e}" for c in logb_model.convert().coef]
    N_model_str = [f"{c:0.4e}" for c in N_model.convert().coef]
    
    
    print("Coefficients for parameter fit:")
    print("logb_param: " + ", ".join(logb_param_str))
    print("N_param: " + ", ".join(N_param_str))
    print("")
    print("Coefficients for model fit:")
    print("logb_model: " + ", ".join(logb_model_str))
    print("N_model: " + ", ".join(N_model_str))
    
    

if __name__ == "__main__":
    # order of polynomial fit for model parameters logb and N
    fit_order = 2
    
    # csv file containing digitized diagram data
    file_path = "ttt_100cr6_kaymak.csv"

    # define the percentage of bainite for each line in the diagram, 
    # in descending order
    line_bainite_percentages = [0.99, 0.01]
    
    # create a pd.DataFrame from the csv
    df = read_csv_file(file_path, line_bainite_percentages)
   
    # note: all times are converted from ms to s for the plots
    
    # reference values from the digitized diagram
    t_ref_1 = df[df["beta"]==0.01]["Time"]*1e-3
    t_ref_99 = df[df["beta"]==0.99]["Time"]*1e-3

    # isothermal fit gives parameters for each unique temperature
    T_i, logb_i, N_i = logb_N_isothermal_fit(df)
    
    # evaluate the parameter fit based on the isothermal fit
    logb_param, N_param = logb_N_parameter_fit(T_i, logb_i, N_i, fit_order)
    
    # evaluate the model fit based on the original data, and based on the
    # parameter fit as the initial guess
    p0 = np.hstack((logb_param.coef, N_param.coef))
    logb_model, N_model = logb_N_model_fit(df, p0)
    
    # text output for polynomial coefficients
    print_coefficients(logb_param, N_param, logb_model, N_model)
    
    # resulting times for 1% and 99% bainite transformation (isothermal fit)
    t_i_1 = jmak_t(0.01, np.exp(logb_i), N_i)*1e-3
    t_i_99 = jmak_t(0.99, np.exp(logb_i), N_i)*1e-3
    
    # choose temperature values for plots of parameter fit
    T_min, T_max = T_i.min(), T_i.max()
    T_param = np.linspace(T_min - 0.1*(T_max - T_min),
                          T_max + 0.1*(T_max - T_min))
    
    # resulting times for 1% and 99% bainite transformation (parameter fit)
    t_param_1 = jmak_t(0.01, np.exp(logb_param(T_param)), N_param(T_param))*1e-3
    t_param_99 = jmak_t(0.99, np.exp(logb_param(T_param)), N_param(T_param))*1e-3

    # choose temperature values for plots of model fit
    T_model = T_param
    
    # resulting times for 1% and 99% bainite transformation (model fit)
    t_model_1 = jmak_t(0.01, np.exp(logb_model(T_model)), N_model(T_model))*1e-3
    t_model_99 = jmak_t(0.99, np.exp(logb_model(T_model)), N_model(T_model))*1e-3

    # plots
    plt.close('all')
    
    # caution: logb_param etc are functions, not arrays, and must be evaluated
    # for the temperatures that will be plotted
    plot_parameters(T_i, logb_i, N_i,
                    T_param, logb_param(T_param), N_param(T_param),
                    T_model, logb_model(T_model), N_model(T_model))
        
    plot_ttt_diagram_all_models(T_i, t_ref_1, t_ref_99, t_i_1, t_i_99, 
                                T_param, t_param_1, t_param_99,
                                T_model, t_model_1, t_model_99)
    
    plot_ttt_diagram_model_fit(T_i, t_ref_1, t_ref_99, t_i_1, t_i_99, 
                               T_model, t_model_1, t_model_99)

