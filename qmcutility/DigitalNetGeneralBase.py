import numpy as np
class DigitalNetGeneralBase:
    def __init__(self, matrices, base):
        """
        self.matrices: list of generating matrices
        self.base: base of the digital net
        self.m: precision of the net
        self.n: log_b(number_of_points)
        self.coeff: coefficient to convert vectors to floats
        """
        self.matrices = matrices[::]
        self.m, self.n = matrices[0].shape
        for M in matrices:
            assert M.shape == (self.m, self.n)
        self.base = base
        self.coeff = np.power(1/self.base, np.arange(1,self.m+1))
    
    def __repr__(self):
        s = f"DigitalNetGeneralBase(dim={len(self.matrices)}, base={self.base})\n"
        s += "Matrices:\n"
        for M in self.matrices:
            s += str(M) + "\n"
        return s

    def generate_all_points(self):
        """
        Generate all points in the unit cube
        """
        points = np.zeros((self.base**self.n,len(self.matrices)))
        for i in range(self.base**self.n):
            vec = self._int_to_vec(i, self.n, self.base)
            for j,Mj in enumerate(self.matrices):
                points[i,j] = (self.coeff@(Mj@vec%self.base)).item()
        return points

    def generate_all_digitally_shifted_points(self,shifts):
        """
        Generate all digitally shifted points in the unit cube
        """
        points = np.zeros((self.base**self.n,len(self.matrices)))

        for i in range(self.base ** self.n):
            vec = self._int_to_vec(i, self.n, self.base)
            for j,Mj in enumerate(self.matrices):
                points[i,j] = (self.coeff@((Mj@vec+shifts[j])%self.base)).item()
        return points

    def kth_vectors(self, k: int):
        vec = self._int_to_vec(k, self.n, self.base)
        return [(Mj@vec%self.base) for j,Mj in enumerate(self.matrices)]

    def kth_vectors_transposed(self, k: int): # transposed for visibility
        return [M.T for M in self.kth_vectors(k)]

    def kth_points(self, k: int):
        return np.array([self.coeff@vec for vec in self.kth_vectors(k)])

    @staticmethod
    def polynomial_lattice(denom, numers, base, n=None, precision=None):
        """
        construct (high order) polynomial lattice
        Args:
          denom: denominator polynomial (as a list of coeff, ascending order)
          numers: list of numerator polynomials (as a list of coeff, ascending order)
          base: base of the lattice
          n: log_b(number of points) (defalt: deg(denom))
          precision: precision of the points (defalt: deg(denom))
        Returns:
          DigitalNetGeneralBase object
        """
        assert denom[-1] == 1 # assume that denom is monic
        if precision == None:
            precision = len(denom)-1
        if n == None:
            n = len(denom)-1

        """
        if q_i(x)/p(x) = (poly) + u_1 x^{-1} + u_2 x^{-2} + ... + u_{m+n-1}x^{m+n-1},
        then the generating matrix C is Hankel with C[i][j] = u_{i+j},
        [[u_1 u_2 u_3 ...]
         [u_2 u_3 u_4 ...]
          ...
        ]
        """
        denom = np.polynomial.Polynomial(denom)
        L = n + precision - 1
        matrices = []
        for numer in numers:
            matrix = np.zeros((precision,n))
            # to avoid padding, concatinate denom
            quotient = (np.polynomial.Polynomial([0]*L + numer)//denom).coef
            result = np.round((np.pad(quotient,(0,L)))[:L][::-1]).astype(np.int64)%base
            for i in range(precision):
                matrix[i,:] = result[i:i+n]
            matrices.append(matrix)

        return DigitalNetGeneralBase(matrices,base)
    
    @staticmethod
    def Faure_sequence(n: int, base: int):
        """
        construct base-dimensional Faure digital net of n*n matrices
        Args:
          n: log_b(number of points)
          base: base of the lattice
        Returns:
          DigitalNetGeneralBase object
        """
        P = DigitalNetGeneralBase.Pascal_upper_matrix(n,base)
        matrices = [np.eye(n, dtype='int'), P]
        for i in range(base-2):
            matrices.append(matrices[-1]@P%base)
        return DigitalNetGeneralBase(matrices,base)
    
    @staticmethod
    def _int_to_vec(N: int, m: int, b: int) -> np.ndarray:
        """
        Convert an integer N to the m-digit vector (np.array) in base b
        """
        res = [0]*m
        for i in range(m):
            N,res[i] = divmod(N,b)
        return np.array(res).reshape(-1, 1)
    
    @staticmethod
    def Pascal_lower_matrix(n: int, base=None):
        P = np.zeros((n,n), dtype='int')
        P[:,0] = 1
        for i in range(1,n):
            for j in range(1,i+1):
                P[i][j] = P[i-1][j-1] + P[i-1][j]
        if base is not None:
            P %= base
        return P

    @staticmethod
    def Pascal_upper_matrix(n: int, base=None):
        return DigitalNetGeneralBase.Pascal_lower_matrix(n,base).T

    @staticmethod
    def identity_matrix(n: int):
        return np.eye(n, dtype='int')
    
    @staticmethod
    def anti_diagonal_identity_matrix(n: int):
        return np.zeros((n,n), dtype='int')[::-1]

    @staticmethod
    def LP_matrix(n: int):
        """
        This function creates an n x n lower triangular matrix where all elements
        above the anti-diagonal, including the diagonal, are set to 1.
        This matrix is used for Larcher-Pillichshammer (0,m,2)-net.
        """
        return (np.tri(n, dtype=int).T)[::-1]

    @staticmethod
    def lower_ones_matrix(n: int):
        """
        This function creates an n x n lower triangular matrix where all elements
        in the lower triangle, including the diagonal, are set to 1
        """
        return np.tri(n, dtype=int)

    @staticmethod
    def random_lower_matrix(n: int, base: int):
        M = np.random.randint(0, base, (n, n))
        for i in range(m):
            M[i,i] = 1
        return np.tril(M)

