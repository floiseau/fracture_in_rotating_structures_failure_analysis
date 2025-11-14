# /// script
# [tool.marimo.display]
# theme = "dark"
# ///

import marimo

__generated_with = "0.16.5"
app = marimo.App(
    app_title="Fracture in rotating structure - Failure prediction",
)

with app.setup:
    # Initialization code that runs before all other cells
    import marimo as mo
    from sympy import symbols, exp, diff, cos, sin, roots, N, Rational, im, pi
    from sympy import latex, sympify


@app.cell(hide_code=True)
def _(K1c, K1e, K2c, K2e, K3c, K3e):
    mo.md(
        rf"""
    # Fracture in rotating structures : Analysis with dimensionless SIFs

    This notebook is dedicated to the analysis of the rotating structure experimental setup.
    It employs the analytical approximation of the SIFs :

    \begin{{align}}
        K_{{I}}^\omega &= {latex(K1c.simplify())}, \\
        K_{{II}}^\omega &= {latex(K2c.simplify())}, \\
        K_{{III}}^\omega &= {latex(K3c.simplify())}, \\
        K_{{I}}^e &= {latex(K1e.simplify())}, \\
        K_{{II}}^e &= {latex(K2e.simplify())}, \\
        K_{{III}}^e &= {latex(K3e.simplify())}.
    \end{{align}}

    ## Determination of failure time

    First we define the different variables in our system.
    """
    )
    return


@app.cell
def _():
    # Define time and its assumptions
    t = symbols("t", real=True, positive=True)

    # Declare the variables with LaTeX names
    Lc_L = symbols(r"\tilde{L}_c", real=True, positive=True)
    a_W = symbols(r"\tilde{a}", real=True, positive=True)
    rho, W, omega, domega, alpha, E, nu, G_c = symbols(
        r"rho W Omega \dot{\Omega} alpha E nu G_c", real=True, positive=True
    )

    # Define derived quantities
    mu = E / (2 * (1 - nu))
    return E, G_c, Lc_L, W, a_W, alpha, domega, mu, nu, omega, rho, t


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    Then, we define the dimensionless SIFs and recombine them into the full SIFs.
    Using the Irwin formula, the energy release rate $G$ is then determined.
    """
    )
    return


@app.cell
def _(E, Lc_L, W, a_W, alpha, domega, mu, nu, omega, rho):
    # Define the dimensionless SIFs (all fractions as rationals)
    K1c_adim = (
        151
        + Rational(171, 10) * exp(a_W) ** 6
        - Rational(112, 10) * exp(a_W) ** 4 * exp(Lc_L) ** 2
    )
    K2c_adim = 0
    K3c_adim = 0
    K1e_adim = (
        1710
        - 3510 * Lc_L * exp(a_W) ** 3
        + 1210 * exp(a_W) ** 4
        + Rational(838, 10) * exp(Lc_L) ** 4
    )
    K2e_adim = -Rational(557, 100) - 127 * a_W**2 - 394 * a_W * (1 - Lc_L)
    K3e_adim = (
        Rational(528, 10)
        + 238 * a_W * (1 - Lc_L)
        + Rational(976, 10) * a_W**3
        - Rational(669, 10) * Lc_L**3
    )

    # Define the SIFs associated with each load
    K1c = rho * omega**2 * W ** (Rational(5, 2)) * K1c_adim
    K2c = rho * omega**2 * W ** (Rational(5, 2)) * K2c_adim
    K3c = rho * omega**2 * W ** (Rational(5, 2)) * K3c_adim
    K1e = rho * domega * W ** (Rational(5, 2)) * K1e_adim * cos(alpha)
    K2e = rho * domega * W ** (Rational(5, 2)) * K2e_adim * cos(alpha)
    K3e = rho * domega * W ** (Rational(5, 2)) * K3e_adim * sin(alpha)

    # Combine the SIFs
    K1 = K1c + K1e
    K2 = K2c + K2e
    K3 = K3c + K3e

    # Computation of the energy release rate from Irwin formula
    Ep = E / (1 - nu**2)
    G = (1 / Ep) * (K1**2 + K2**2) + K3**2 / (2 * mu)
    return G, K1, K1c, K1e, K2, K2c, K2e, K3, K3c, K3e


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""Finally, the failure time is given by Griffith criterion $G=G_c$. As the final expression of failure time depends on a lot of paramaters, and is not explicitly required, we directly substitute the parameters by there numerical values."""
    )
    return


@app.cell(hide_code=True)
def _():
    subs_values = mo.md(
        r"""
        The material parameters are
        - Density: $\rho=${rho} kg/m$^3$,
        - Young modulus: $E=${E} Pa,
        - Poisson ratio: $\nu=${nu},
        - Fracture toughness: $G_c=${G_c} J/m$^2$.

        The geometric parameters are 
        - Plate width: $W=$ {W} m,
        - Crack length: $\frac{{a}}{{W}}=$ {a_W},
        - Crack position: $\frac{{L_c}}{{L}}=$ {Lc_L}.

        The load is
        - $\Omega(t) =$ {omega} rad/s,
        - $\alpha =$ {alpha} deg.
        """
    ).batch(
        # Material
        rho=mo.ui.number(1200),
        E=mo.ui.number(3_200_000),
        nu=mo.ui.number(0.4),
        G_c=mo.ui.number(220),
        # Geometry
        W=mo.ui.number(40e-3),
        a_W=mo.ui.number(1 / 2),
        Lc_L=mo.ui.number(1 / 2),
        # Load
        omega=mo.ui.text("10 * t"),
        alpha=mo.ui.number(0),
    )
    subs_values
    return (subs_values,)


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    In the following cell, we compute the failure time for a specific configuration.
    Then, the load properties at this time are shown (SIFs, mode mixity, etc).
    The user must define a load (either in velocity or in acceleration), along with all the geometrical, and mechanical properties.
    """
    )
    return


@app.cell
def _(
    E,
    G,
    G_c,
    K1,
    K2,
    K3,
    Lc_L,
    W,
    a_W,
    alpha,
    domega,
    nu,
    omega,
    rho,
    subs_values,
    t,
):
    # subs = {
    #     # Material
    #     rho: 1200,
    #     E: 3200000000,
    #     nu: Rational(4, 10),
    #     G_c: 2700,
    #     # Geometry
    #     a_W: Rational(1, 2),
    #     Lc_L: Rational(1, 2),
    #     W: Rational(40, 1000),
    #     # Load
    #     omega: 900 * t,
    #     domega: diff(900 * t, t),
    #     alpha: 0,
    # }

    # Gather the parameters
    str_to_var = {
        "rho": rho,
        "E": E,
        "nu": nu,
        "G_c": G_c,
        "a_W": a_W,
        "Lc_L": Lc_L,
        "W": W,
        "omega": omega,
        "alpha": alpha,
    }
    subs = {str_to_var[key]: val for key, val in subs_values.value.items()}
    subs[omega] = sympify(subs[omega], locals={"t": t})
    subs[domega] = diff(subs[omega], t)
    subs[alpha] *= 180 / pi

    # Determine the failure time
    failure_time_equation = (G - G_c).subs(subs).simplify()
    print("The failure time equation (based on Griffith criterion) is:")
    print(latex(failure_time_equation.simplify()))

    sols = roots(failure_time_equation, t)
    t_sols = [t_sol.subs(subs) for t_sol in sols]
    print("\nThe possible failure times are:")
    for i, t_sol in enumerate(t_sols):
        print(f"t_{i + 1} ≈ {N(t_sol)} seconds")

    # # Extract the failure time
    try:
        valid_ts = [N(t_sol) for t_sol in t_sols if im(N(t_sol)) == 0 and N(t_sol) >= 0]
        t_failure = min(valid_ts)
    except:
        # If all the values are complex or negative, it means that the acceleration load is sufficient to break the beam at the very begining of the load.
        print("The initial load is sufficient for failure => t_failure = 0.")
        t_failure = 0
    subs[t] = t_failure

    # Determine the load at failure
    omega_failure = omega.subs(subs)
    domega_failure = domega.subs(subs)

    # Determine the energy release rate at failure
    G_failure = G.subs(subs)

    # Determine the associated SIFs
    K1_failure = K1.subs(subs)
    K2_failure = K2.subs(subs)
    K3_failure = K3.subs(subs)

    # TODO: Calculate the initial crack propagation angle.
    print("\nThen, we could predict the initial crack propagation angle.")
    return (
        G_failure,
        K1_failure,
        K2_failure,
        K3_failure,
        domega_failure,
        omega_failure,
        subs,
        t_failure,
    )


@app.cell(hide_code=True)
def _(
    G_failure,
    K1_failure,
    K2_failure,
    K3_failure,
    domega_failure,
    omega_failure,
    t_failure,
):
    mo.md(
        rf"""
    Based on the previous calculations, we obtain:

    - the failure time :
        - $t_f = {t_failure:.3g}$ s,
    - the load at failure :
        - $\Omega = {float(omega_failure):.3g}$ rad/s,
        - $\dot{{\Omega}} = {float(domega_failure):.3g}$ rad/s$^2$,
    - the energy release rate at failure :
        - $G = {G_failure:.3g}$ J,
    - the SIFs at failure :
        - $K_{{I}} = {float(K1_failure):.3g}$ Pa.m$^{{1/2}}$
        - $K_{{II}} = {float(K2_failure):.3g}$ Pa.m$^{{1/2}}$
        - $K_{{III}} = {float(K3_failure):.3G}$ Pa.m$^{{1/2}}$
    - the mode mixities :
        - $M_{{I}}^{{II}} = ...$
        - $M_{{I}}^{{III}} = ...$
    - the initial crack propagation angle :
        - $\varphi = ...$
    """
    )
    return


@app.cell(hide_code=True)
def _(K1, K2, subs, t, t_failure):
    import matplotlib.pyplot as plt
    import numpy as np

    new_subs = {}
    new_subs.update(subs)
    if t in new_subs:
        del new_subs[t]

    t_max = t_failure if t_failure > 0 else 1
    time_values = np.linspace(0, t_max, 100)
    K1_values = [K1.subs(new_subs).subs({t: time}) for time in time_values]
    K2_values = [K2.subs(new_subs).subs({t: time}) for time in time_values]
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(time_values, K1_values, label="$K_{I}$")
    plt.plot(time_values, K2_values, label="$K_{II}$")
    plt.axvline(t_failure, c="red", label="$t_c$")
    plt.xlabel("Time (s)")
    plt.ylabel("$K_p$ (Pa.m$^{1/2}$)")
    plt.title("SIF as a function of time")
    plt.grid(True)
    plt.legend()
    plt.gca()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
