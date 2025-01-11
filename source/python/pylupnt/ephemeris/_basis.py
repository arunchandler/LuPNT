import numpy as np
from abc import ABC, abstractmethod 

class Basis(ABC):
    """
    Abstract class for basis functions.
    """

    @abstractmethod
    def generate_pol(self, n, t):
        """
        Generate the polynomial basis for the given degree n and time vector t.

        Args:
            n (int): polynomial order
            t (float): time value

        Returns:
            loat: coefficient of order n at time t

        """
        pass

    @abstractmethod
    def generate_der_pol(self, n, t, t_interval):
        """
        Generate the derivative of the polynomial basis for the given degree n and time vector t.

        Args:
            n (int): polynomial order
            t (float): time value

        Returns:
            float: derivative of coefficient of order n at time t
        """
        pass

    def get_f_approx(self, a, n, t_array):
        """
        Get the approximation of the function f(t) using the basis functions.

        Args:
            a (np.array): coefficients of the basis functions
            n (int): polynomial order
            t_array (np.array): time vector

        Returns:
            np.array: approximation of f(t)
        """
        f_approx = np.zeros(len(t_array))

        # TODO: let's make this a lot faster?

        for t in range(len(t_array)):
            for n_val in range(n + 1):
                f_approx[t] += a[n_val] * self.generate_pol(n_val, t_array[t])

        return f_approx

    def get_f_approx_der(self, a, n, t_array, t_interval):
        """
        Get the approximation of the derivative of the function f(t) using the basis functions.

        Args:
            a (np.array): coefficients of the basis functions
            n (int): polynomial order
            t_array (np.array): time vector

        Returns:
            np.array: approximation of f'(t)
        """
        f_approx = np.zeros(len(t_array))
        for t in range(len(t_array)):
            for n_val in range(n + 1):
                f_approx[t] += a[n_val] * self.generate_der_pol(
                    n_val, t_array[t], t_interval
                )

        return f_approx


############################################
# Chebyshev basis functions
############################################
class Chebyshev(Basis):
    """
    Chebyshev basis functions.
    """

    def generate_pol(self, n, t):
        """Generates Chebyshev coefficients recursively

        Args:
            n (int): polynomial order
            t (float): time value

        Returns:
            float: coefficient of order n at time t
        """
        if n == 0:
            return 1
        if n == 1:
            return t
        if n >= 2:
            return 2 * t * (self.generate_pol(n - 1, t)) - (
                self.generate_pol(n - 2, t)
            )

    def generate_der_pol(self, n, t, t_interval):
        """Generates derivative of Chebyshev coefficients recursively

        Args:
            n (int): polynomial order
            t (float): time value
            t_interval (float): total time interval

        Returns:
            float: derivative of coefficient of order n at time t
        """
        if n == 0:
            return 0
        if n == 1:
            return 1 * (2 / t_interval)
        if n >= 2:
            # only place 2/t_interval at the generate_pol, so that as a total,
            # the derivative is multiplied by 2/t_interval
            return (
                2 * t * (self.generate_der_pol(n - 1, t, t_interval))
                + 2 * (self.generate_pol(n - 1, t)) * (2 / t_interval)
                - (self.generate_der_pol(n - 2, t, t_interval))
            )

    def get_vel_coeffs(self, n, a):
        """Calculate derivative approximation coefficients based on function approximation coefficients

        Args:
            n (int): polynomial order
            a (list): list of approximation coefficients

        Returns:
            float: derivative coefficients not using derivative of Cheybeshev polynomials
        """
        if n == 0:
            return a[1] + self.get_vel_coeffs(2, a) / 2

        if n >= 1 and n <= len(a) - 3:
            return 2 * (n + 1) * (a[n + 1]) - self.get_vel_coeffs(n + 2, a)

        if n == len(a) - 2:
            return 2 * (len(a) - 1) * a[-1]

        if n == len(a) - 1:  # zero indexed
            return 0
        

############################################
# Polynomial basis functions
############################################
class Polynomial(Basis):
    """
    Polynomial basis functions.
    """

    def generate_pol(self, n, t):
        """Generates polynomial coefficients

        Args:
            n (int): polynomial order
            t (float): time value

        Returns:
            float: coefficient of order n at time t
        """
        return t ** n
    
    def generate_der_pol(self, k, t, t_interval):
        """
        Generates k-th order polynomial derivative evaluation at time t

        Args:
            k (int): polynomial order
            t (float): time value
            t_interval (float): total time interval

        Returns:
            float: derivative component of order k at time t
        """
        if k == 0:
            return 0
        else:
            return k * (t ** (k - 1)) * (2 / t_interval)


############################################
# Legrendre basis functions
############################################
class Legendre(Basis):
    """
    Legendre basis functions.
    """

    def generate_pol(self, n, t):
        """Generates Legendre coefficients

        Args:
            n (int): polynomial order
            t (float): normalized time value [-1, 1]

        Returns:
            float: coefficient of order n at time t
        """
        if n == 0:
            return 1
        if n == 1:
            return t
        if n >= 2:
            return ((2 * n - 1) * t * self.generate_pol(n - 1, t) 
                    - (n - 1) * self.generate_pol(n - 2, t)) / n
        
    
    def generate_der_pol(self, n, t, t_interval):
        """
        Generates derivative of Legendre coefficients

        Args:
            n (int): polynomial order
            t (float): normalized time value [-1, 1]
            t_interval (float): total time interval

        Returns:
            float: derivative of coefficient of order n at time t
        """
        if n == 0:
            return 0
        if n == 1:
            return 1.0 * (2 / t_interval)
        if n >= 2:
            # only place 2/t_interval at the generate_pol, so that as a total,
            # the derivative is multiplied by 2/t_interval
            return ((2 * n - 1) * self.generate_pol(n - 1, t) * (2 / t_interval) 
                      + (2 * n - 1) * t * self.generate_der_pol(n - 1, t, t_interval) 
                      - (n - 1) * self.generate_der_pol(n - 2, t, t_interval)) / n


############################################
# Fourier basis functions
############################################

class Fourier(Basis):
    """
    Fourier basis functions + Linear term
    """

    def generate_pol(self, n, t):
        """Generates Fourier coefficients

        Args:
            n (int): polynomial order
            t (float): time value

        Returns:
            float: coefficient of order n at time t
        """
        if n == 0:
            return 1
        elif n == 1:
            return t
        else:
            k = n // 2
            if n % 2 == 0:
                # cos(k * pi * t)
                return np.cos(k * np.pi * t)
            else:
                # sin(k * pi * t)
                return np.sin(k * np.pi * t)
        
    
    def generate_der_pol(self, n, t, t_interval):
        """
        Generates derivative of Fourier coefficients

        Args:
            n (int): polynomial order
            t (float): time value
            t_interval (float): total time interval

        Returns:
            float: derivative of coefficient of order n at time t
        """
        if n == 0:
            return 0
        elif n == 1:
            return 1.0 * (2 / t_interval)
        else:
            k = n // 2
            if n % 2 == 0:
                # cos(k * pi * t) => derivative is -k*pi sin(k*pi*t)
                return -k * np.pi * np.sin(k * np.pi * t) * (2/t_interval)
            else:
                # sin(k * pi * t) => derivative is k*pi cos(k*pi*t)
                return  k * np.pi * np.cos(k * np.pi * t) * (2/t_interval)