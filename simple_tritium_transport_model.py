import pint
import numpy as np
from scipy.optimize import fsolve

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

        self.L_wall = 1 * ureg.cm  # TODO: VERIFY THIS

        self.neutron_rate = 3e8 * ureg.neutron * ureg.s**-1

        self.c_old = 0 * ureg.particle * ureg.m**-3
        self.dt = 1 * ureg.h
        self.exposure_time = 300 * ureg.h
        self.TBR = TBR

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
        return perimeter_wall * self.height

    def source(self, t):
        if t < self.exposure_time:
            return self.TBR * self.neutron_rate
        else:
            return 0 * self.TBR * self.neutron_rate

    def Q_wall(self, c_salt):
        k_wall = 1.9e-8 * ureg.m * ureg.s**-1  # from Kumagai
        return self.A_wall * k_wall * c_salt

    def Q_top(self, c_salt):
        k_top = 4.9e-7 * ureg.m * ureg.s**-1  # from Kumagai
        return self.A_top * k_top * c_salt

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


def quantity_to_activity(Q):
    return Q * SPECIFIC_ACT * MOLAR_MASS
