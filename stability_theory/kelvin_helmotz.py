import cmath


class KelvinHelmholtz:
    def __init__(self, rho1, rho2, U1, U2, g=0, sigma=0):
        self.rho1 = rho1
        self.rho2 = rho2
        self.U1 = U1
        self.U2 = U2
        self.g = g
        self.sigma = sigma

    def complex_wave_number(self, k) -> complex:
        """
        Calculate the Kelvin-Helmholtz wave number for a two-layer fluid system.

        Parameters:
        rho1 : float
            Density of the upper fluid layer.
        rho2 : float
            Density of the lower fluid layer.
        U1 : float
            Velocity of the upper fluid layer.
        U2 : float
            Velocity of the lower fluid layer.
        g : float, optional
            Gravitational acceleration (default is 0).
        sigma : float, optional
            Surface tension (default is 0).

        Returns:
        k : float
            Kelvin-Helmholtz wave number.
        """
        rho_p = self.rho1 + self.rho2
        rho_m = self.rho1 - self.rho2

        # take greatest solution of the quadratic eq, which in this particular case corresponds to the positive branch
        c = ((self.rho1*self.U1 + self.rho2*self.U2) / rho_p) + cmath.sqrt(
            (rho_m*self.g)/(rho_p*k)  # gravitational term
            + (self.sigma*k)/rho_p  # surface tension term
            - (self.rho1*self.rho2*(self.U1 - self.U2)**2) /
            (rho_p**2)  # "shear" term
        )
        return c

    def growth_rate(self, k) -> float:
        """
        Calculate the growth rate of the Kelvin-Helmholtz instability for a given wave number.

        Parameters:
        k : float
            Wave number.

        Returns:
        float
            Growth rate of the instability.
        """
        c = self.complex_wave_number(k)
        return k*c.imag

    def frequency(self, k) -> float:
        """
        Calculate the frequency of the Kelvin-Helmholtz instability for a given wave number.

        Parameters:
        k : float
            Wave number.

        Returns:
        float
            Frequency of the instability.
        """
        c = self.complex_wave_number(k)
        return k*c.real


if __name__ == "__main__":
    # Example usage
    k = 1.0
    print("Initial parameters:")
    print("rho1=1.0, rho2=2.0, U1=1.0, U2=0.5, g=9.81, sigma=0.07")
    kh = KelvinHelmholtz(rho1=1.0, rho2=2.0, U1=1.0,
                         U2=0.5, g=9.81, sigma=0.07)
    print("Complex wave number:", kh.complex_wave_number(k))
    print("Growth rate:", kh.growth_rate(k))
    print("Frequency:", kh.frequency(k))

    # test U1 and U2 = 0 (capillary-gravity wave)
    print("\nTesting zero velocity (U1=U2=0):")
    print("rho1=1.0, rho2=2.0, U1=0.0, U2=0.0, g=9.81, sigma=0.07")
    kh_zero_velocity = KelvinHelmholtz(
        rho1=1.0, rho2=2.0, U1=0.0, U2=0.0, g=9.81, sigma=0.07)
    print("Complex wave number (U1=U2=0):",
          kh_zero_velocity.complex_wave_number(k))
    print("Growth rate (U1=U2=0):", kh_zero_velocity.growth_rate(k))
    print("Frequency (U1=U2=0):", kh_zero_velocity.frequency(k))

    # test: absence of gravity and surface tension (pure shear instability)
    print("\nTesting pure shear instability (g=sigma=0):")
    print("rho1=1.0, rho2=2.0, U1=1.0, U2=0.5, g=0.0, sigma=0.0")
    kh_no_gravity_no_surface = KelvinHelmholtz(
        rho1=1.0, rho2=2.0, U1=1.0, U2=0.5, g=0.0, sigma=0.0
    )
    print("Complex wave number (g=sigma=0):",
          kh_no_gravity_no_surface.complex_wave_number(k))
    print("Growth rate (g=sigma=0):", kh_no_gravity_no_surface.growth_rate(k))
    print("Frequency (g=sigma=0):", kh_no_gravity_no_surface.frequency(k))
