import numpy as np
import matplotlib.pyplot as plt
from plotly import graph_objects as go

# Global variables
ALPHA = 100
UNINFECTED_PRODUCTION = 0.272
UNINFECTED_DEATH_RATE = 0.00136
LATENT_DEATH_RATE = UNINFECTED_DEATH_RATE
INFECTED_DEATH_RATE = 0.33
VIRION_INFECTION_RATE = 0.00027
VIRUS_DEATH_RATE = 2
TREATMENT_KILL_RATE = 2

UNINFECTED_TO_LATENT_RATE = 0.2
UNINFECTED_TO_INFECTED_RATE = 1 - UNINFECTED_TO_LATENT_RATE

RUN_TIME = 7500  # days
TIME_STEP = 0.05  # days
NUM_STEPS = int(RUN_TIME / TIME_STEP)
LATENT_BEGINNING = 2000 / TIME_STEP

def delta_u(uninfected, virion):
    return ((UNINFECTED_PRODUCTION - uninfected * virion * VIRION_INFECTION_RATE - UNINFECTED_DEATH_RATE * uninfected) *
            TIME_STEP)


def delta_i(uninfected, virion, latent, infected, latent_rate=0.05):
    return ((UNINFECTED_TO_INFECTED_RATE * uninfected * virion * VIRION_INFECTION_RATE) + (latent_rate * latent)
            - (INFECTED_DEATH_RATE * infected)) * TIME_STEP


def delta_v(virion, infected, treatment=False):
    if not treatment:
        return (ALPHA * infected - VIRUS_DEATH_RATE * virion) * TIME_STEP
    else:
        return (ALPHA * infected - VIRUS_DEATH_RATE * virion - TREATMENT_KILL_RATE * virion) * TIME_STEP


def delta_l(uninfected, virion, latent, latent_rate=0.05):
    return (UNINFECTED_TO_LATENT_RATE * uninfected * virion * VIRION_INFECTION_RATE - LATENT_DEATH_RATE * latent -
            latent_rate * latent) * TIME_STEP


def run_step(virion, uninfected, latent, infected, latent_rate, i, treatment):
    virion[i] = virion[i - 1] + delta_v(virion[i - 1], infected[i - 1], treatment)
    uninfected[i] = uninfected[i - 1] + delta_u(uninfected[i - 1], virion[i - 1])
    latent[i] = latent[i - 1] + delta_l(uninfected[i - 1], virion[i - 1], latent[i - 1],  latent_rate)
    infected[i] = infected[i - 1] + delta_i(uninfected[i - 1], virion[i - 1], latent[i - 1], infected[i - 1],
                                            latent_rate)


def original_model(virion, uninfected, latent, infected):
    # Run twice, once with latent to infected from begging, and once with latent to infected after reaching equilibrium
    for i in range(2):
        for j in range(1, NUM_STEPS):
            if i == 0:  # If latent to infected starts from beginning
                latent_rate = 0.05
                run_step(virion, uninfected, latent, infected, latent_rate, j, False)
            else:  # If latent outbreaks only after reaching equilibrium
                if j < LATENT_BEGINNING:
                    latent_rate = 0
                    run_step(virion, uninfected, latent, infected, latent_rate, j, False)
                else:
                    latent_rate = 0.05
                    run_step(virion, uninfected, latent, infected, latent_rate, j, False)

        # Plot all the equations on same graph with plotly
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=np.arange(0, RUN_TIME, TIME_STEP), y=virion, name="Virion"))
        fig.add_trace(go.Scatter(x=np.arange(0, RUN_TIME, TIME_STEP), y=uninfected, name="Uninfected"))
        fig.add_trace(go.Scatter(x=np.arange(0, RUN_TIME, TIME_STEP), y=latent, name="Latent"))
        fig.add_trace(go.Scatter(x=np.arange(0, RUN_TIME, TIME_STEP), y=infected, name="Infected"))
        if i == 0:  # If latent to infected starts from beginning
            fig.update_layout(title="HIV model with latent to infected from beginning", xaxis_title="Time (days)",
                              yaxis_title="Number of cells")
            fig.write_html(f"../Figs/HIV Model Latent from Beginning.html", full_html=False, include_plotlyjs='cdn')

        else:  # If latent outbreaks only after reaching equilibrium
            # Add vertical line to show when latent to infected starts
            fig.add_shape(
                dict(
                    type='line',
                    x0=2000,
                    x1=2000,
                    y0=0,
                    y1=12500,
                    line=dict(color='black', width=2, dash='dash'),
                )
            )

            # Add text annotation
            fig.add_annotation(
                dict(
                    x=1800,
                    y=12000,
                    text='Latent Outbreak',
                    showarrow=False,
                    font=dict(color='black')
                )
            )
            fig.update_layout(title="HIV model with latent to infected after reaching equilibrium",
                              xaxis_title="Time (days)",
                              yaxis_title="Number of cells")
            fig.write_html(f"../Figs/HIV Model Latent after Equilibrium.html", full_html=False, include_plotlyjs='cdn')

        fig.show()


def treatment_model(virion, uninfected, latent, infected):
    # Run twice, once with just treatment, and once with treatment and higher latent to infected rate
    for i in range(2):
        for j in range(1, NUM_STEPS):
            if i == 0:  # If just treatment
                if j < LATENT_BEGINNING:
                    latent_rate = 0
                else:
                    latent_rate = 0.05
                run_step(virion, uninfected, latent, infected, latent_rate, j, True)
            else:  # If treatment and higher latent to infected rate
                if j < LATENT_BEGINNING:
                    latent_rate = 0
                else:
                    latent_rate = 0.5
                run_step(virion, uninfected, latent, infected, latent_rate, j, True)

        # Plot all the equations on same graph with plotly
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=np.arange(0, RUN_TIME, TIME_STEP), y=virion, name="Virion"))
        fig.add_trace(go.Scatter(x=np.arange(0, RUN_TIME, TIME_STEP), y=uninfected, name="Uninfected"))
        fig.add_trace(go.Scatter(x=np.arange(0, RUN_TIME, TIME_STEP), y=latent, name="Latent"))
        fig.add_trace(go.Scatter(x=np.arange(0, RUN_TIME, TIME_STEP), y=infected, name="Infected"))
        if i == 0:
            fig.update_layout(title="HIV model with treatment", xaxis_title="Time (days)",
                              yaxis_title="Number of cells")
            fig.write_html(f"../Figs/HIV Model Treatment.html", full_html=False, include_plotlyjs='cdn')
        else:
            fig.update_layout(title="HIV model with treatment and higher latent to infected rate",
                              xaxis_title="Time (days)",
                              yaxis_title="Number of cells")
            fig.write_html(f"../Figs/HIV Model Treatment and Latent.html", full_html=False, include_plotlyjs='cdn')
        fig.show()


def init_arrays():
    virion = np.zeros(NUM_STEPS)
    uninfected = np.zeros(NUM_STEPS)
    latent = np.zeros(NUM_STEPS)
    infected = np.zeros(NUM_STEPS)

    virion[0] = 10
    uninfected[0] = 500

    return virion, uninfected, latent, infected


if __name__ == '__main__':
    virion, uninfected, latent, infected = init_arrays()
    original_model(virion, uninfected, latent, infected)

    virion, uninfected, latent, infected = init_arrays()
    treatment_model(virion, uninfected, latent, infected)
