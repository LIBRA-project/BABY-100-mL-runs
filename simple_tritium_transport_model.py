import pint
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import cumulative_trapezoid

ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
ureg.define("neutron = 1 * particle = n")

SPECIFIC_ACT = 3.57e14 * ureg.Bq * ureg.g**-1
MOLAR_MASS = 6.032 / 2 * ureg.g * ureg.mol**-1


class Model:
    def __init__(
        self,
        radius: pint.Quantity = None,
        height: pint.Quantity = None,
        TBR: pint.Quantity = None,
    ) -> None:
        self.radius = radius
        self.height = height

        self.L_wall = 0.06 * ureg.inches

        self.neutron_rate = 3e8 * ureg.neutron * ureg.s**-1

        self.c_old = 0 * ureg.particle * ureg.m**-3
        self.dt = 1 * ureg.h
        self.exposure_time = 300 * ureg.h
        self.number_days = 1 * ureg.day
        self.TBR = TBR

        self.k_wall = 1.9e-8 * ureg.m * ureg.s**-1  # from Kumagai
        self.k_top = 4.9e-7 * ureg.m * ureg.s**-1  # from Kumagai

        self.concentrations = []
        self.times = []

    @property
    def volume(self):
        return self.A_top * self.height

    @property
    def A_top(self):
        return np.pi * self.radius**2

    @property
    def A_wall(self):
        perimeter_wall = 2 * np.pi * (self.radius + self.L_wall)
        return perimeter_wall * self.height + self.A_top

    def source(self, t):
        if t % (24 * ureg.h) < self.exposure_time and t < self.number_days:
            return self.TBR * self.neutron_rate
        else:
            return 0 * self.TBR * self.neutron_rate

    def Q_wall(self, c_salt):
        return self.A_wall * self.k_wall * c_salt

    def Q_top(self, c_salt):
        return self.A_top * self.k_top * c_salt

    def equation(self, c, t):
        c *= ureg.particle * ureg.m**-3

        lhs = self.volume * (c - self.c_old) / self.dt
        rhs = self.source(t) - self.Q_wall(c) - self.Q_top(c)
        return lhs - rhs

    def run(self, t_final):
        t = 0 * ureg.s
        while t < t_final:
            t += self.dt
            c_new = (
                fsolve(self.equation, x0=self.c_old, args=(t,))
                * ureg.particle
                * ureg.m**-3
            )
            self.c_old = c_new

            self.times.append(t)
            self.concentrations.append(c_new)
        self.concentrations = ureg.Quantity.from_list(self.concentrations)
        self.times = ureg.Quantity.from_list(self.times)

    def reset(self):
        self.c_old = 0 * ureg.particle * ureg.m**-3
        self.concentrations = []
        self.times = []

    def integrated_release_top(self):
        top_release = self.Q_top(self.concentrations)
        integrated_top = cumulative_trapezoid(
            top_release.to(ureg.particle * ureg.h**-1),
            self.times.to(ureg.h),
            initial=0,
        )
        integrated_top *= ureg.particle  # attach units
        return integrated_top

    def integrated_release_wall(self):
        wall_release = self.Q_wall(self.concentrations)
        integrated_wall = cumulative_trapezoid(
            wall_release.to(ureg.particle * ureg.h**-1),
            self.times.to(ureg.h),
            initial=0,
        )
        integrated_wall *= ureg.particle  # attach units
        return integrated_wall


def quantity_to_activity(Q):
    return Q * SPECIFIC_ACT * MOLAR_MASS


def activity_to_quantity(A):
    return A / (SPECIFIC_ACT * MOLAR_MASS)