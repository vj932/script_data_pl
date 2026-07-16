
#### Data Provided

- Survival curve of both antibiotics over time (in CFU with respect to time)

- Growth curves of transformed bacteria (in $OD_{600}$)

- Final proportion of plasmid-containing colonies to plasmid-free colonies upon plating

  

#### Assumptions made

- Sister cells of bacteria which no longer contain their plasmid (F) still have the same survival probability as parent

- Plasmid loss rate is a constant and does not vary according to other factors like time since inoculation or relative population density

- PFB do not replicate under antibiotic exposure

- The $OD_{600}$ curve provided can be directly used to find doubling time

- Plasmid-containing bacteria (C) do not die over the course of the experiment

- Plasmid loss happens at division

  
  
  
  

#### Derivation of solution:

Where $\beta$ is the growth rate, $C_0, F_0$ are the initial proportions of C and F, L is the plasmid loss rate with respect to time, $t$ is the variable time, $T$ is the length of the survival curve(determined by maximum bacteria survival time for a treatment), and p(t) is the probability of a plasmid-free surviving at $T - t$ hours away from plating.

  

The change in C with respect to time is proportionate to its growth minus its loss, both of which are proportionate to C itself, which motivates the differential equation:$\\$

$$\frac{dC}{dt} = C*\beta - C * L = C(\beta - L) = C*(\beta - L)$$
Solving for C given the DE above gives the following: $\\$

$$C(t) = C_0 * e^{\left((\beta - L) * t \right)}$$
Which also gives us $\beta$ looking at just the growth portion of the equation. Since division is doubling, and we have the doubling time $d$ from the data:$\\$

$$

\begin{aligned}

P &= e^{\beta * t}\\

2 &= e^{\beta * d}\\

\frac{\ln 2}{d} &= \beta

\end{aligned}

$$
The change in the plasmid free proportion of bacteria is proportionate to the loss rate, probability of survival at time of loss, and the solution for C which we derived above. $\\$

$$\frac{dF}{dt} = C(t) * L * p(t) = e^{\left((\beta - L) * t \right)} * L * p(t)$$
Integration and substitution of $C(t)$ gives the following expression: $\\$

$$F = L * \int_{0}^{T}e^{\left((\beta - L) * t \right)} * p(t)dt $$
And since the final piece of data gives us the proportion of $R = \frac{C}{F}$, and since only the proportion of $\frac{C}{F}$ matter, we can set the initials $C_0 = 1, F_0 = 0$ as C represents 100% of the population at t = 0 which appears at plating. As a result they do not appear in the final equation: $\\$

$$R = \frac{e^{(\beta - L) * t}}{L * \int_{0}^{T}e^{(\beta - L) * t} * p(t)dt}$$
And as L cannot be isolated in the expression above, it has to be computed numerically. Since $L$ is initially computed with respect to time, and loss per division, $L_{div}$ is more biologically interpretable, we multiply $L$ by the doubling time. $$L_{div} = L * d$$