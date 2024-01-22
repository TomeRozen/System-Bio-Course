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
LATENT_TO_INFECTED_RATE = 0.05
VIRUS_DEATH_RATE = 2  # TODO: Make sure it makes sense when deciding on time step

UNINFECTED_TO_LATENT_RATE = 0.2
UNINFECTED_TO_INFECTED_RATE = 1 - UNINFECTED_TO_LATENT_RATE

RUN_TIME = 3500  # days
TIME_STEP = 0.05  # days
NUM_STEPS = int(RUN_TIME / TIME_STEP)
LATENT_BEGINNING = 2000 / TIME_STEP

def delta_u(uninfected, virion):
    return ((UNINFECTED_PRODUCTION - uninfected * virion * VIRION_INFECTION_RATE - UNINFECTED_DEATH_RATE * uninfected) *
            TIME_STEP)


def delta_i(uninfected, virion, latent, infected):
    return ((UNINFECTED_TO_INFECTED_RATE * uninfected * virion * VIRION_INFECTION_RATE) + (LATENT_TO_INFECTED_RATE * latent)
            - (INFECTED_DEATH_RATE * infected)) * TIME_STEP


def delta_v(virion, infected):
    return (ALPHA * infected - VIRUS_DEATH_RATE * virion) * TIME_STEP


def delta_l(uninfected, virion, latent):
    return (UNINFECTED_TO_LATENT_RATE * uninfected * virion * VIRION_INFECTION_RATE - LATENT_DEATH_RATE * latent -
            LATENT_TO_INFECTED_RATE * latent) * TIME_STEP


def run_step(virion, uninfected, latent, infected, i):
    virion[i] = virion[i - 1] + delta_v(virion[i - 1], infected[i - 1])
    uninfected[i] = uninfected[i - 1] + delta_u(uninfected[i - 1], virion[i - 1])
    latent[i] = latent[i - 1] + delta_l(uninfected[i - 1], virion[i - 1], latent[i - 1])
    infected[i] = infected[i - 1] + delta_i(uninfected[i - 1], virion[i - 1], latent[i - 1], infected[i - 1])



if __name__ == '__main__':

    # Run twice, once with latent to infected from begging, and once with latent to infected after reaching equilibrium
    for i in range(2):
        virion = np.zeros(NUM_STEPS)
        uninfected = np.zeros(NUM_STEPS)
        latent = np.zeros(NUM_STEPS)
        infected = np.zeros(NUM_STEPS)

        virion[0] = 10
        uninfected[0] = 500

        for j in range(1, NUM_STEPS):
            if i == 0:
                run_step(virion, uninfected, latent, infected, j)
            else:
                if j < LATENT_BEGINNING:
                    LATENT_TO_INFECTED_RATE = 0
                else:
                    LATENT_TO_INFECTED_RATE = 0.05
                run_step(virion, uninfected, latent, infected, j)

        # Plot all the equations on same graph with plotly
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=np.arange(0, RUN_TIME, TIME_STEP), y=virion, name="Virion"))
        fig.add_trace(go.Scatter(x=np.arange(0, RUN_TIME, TIME_STEP), y=uninfected, name="Uninfected"))
        fig.add_trace(go.Scatter(x=np.arange(0, RUN_TIME, TIME_STEP), y=latent, name="Latent"))
        fig.add_trace(go.Scatter(x=np.arange(0, RUN_TIME, TIME_STEP), y=infected, name="Infected"))
        if i == 0:
            fig.update_layout(title="HIV model with latent to infected from beginning", xaxis_title="Time (days)",
                              yaxis_title="Number of cells")
            fig.write_html(f"../Figs/HIV Model Latent from Beginning.html")
        else:
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
            fig.update_layout(title="HIV model with latent to infected after reaching equilibrium", xaxis_title="Time (days)",
                              yaxis_title="Number of cells")
            fig.write_html(f"../Figs/HIV Model Latent after Equilibrium.html")

        fig.show()
