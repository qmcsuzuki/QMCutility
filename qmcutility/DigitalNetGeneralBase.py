class DigitalNetGeneralBase:
    def __init__(self, matrices, base):
        self.matrices = matrices[::]
        self.dim = len(matrices)
        self.m, self.n = matrices[0].shape # m: precision, n: log_b(number_of_points)
        for M in matrices:
            assert M.shape == (self.m, self.n)
        self.base = base
    
    def generate_all_points(self):
        """
        Generate all points in the unit cube
        """
        ratio = np.power(1/self.base, np.arange(1,self.m+1))

        points = np.zeros((self.base**self.n,len(self.matrices)))
        for i in range(self.base**self.n):
            vec = self._int_to_vec(i, self.n, self.base)
            for j,Mj in enumerate(self.matrices):
                points[i,j] = (ratio@(Mj@vec%self.base)).item()
        return points

    def generate_all_digitally_shifted_points(self,shifts):
        """
        Generate all digitally shifted points in the unit cube
        """
        ratio = np.power(1/self.base, np.arange(1,self.m+1))

        points = np.zeros((self.base**self.n,len(self.matrices)))

        for i in range(self.base ** self.n):
            vec = self._int_to_vec(i, self.n, self.base)
            for j,Mj in enumerate(self.matrices):
                points[i,j] = (ratio@((Mj@vec+shifts[j])%self.base)).item()
        return points

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
    def _int_to_vec(N: int, m: int, b: int) -> np.ndarray:
        """
        Convert an integer N to the m-digit vector (np.array) in base b
        """
        res = [0]*m
        for i in range(m):
            N,res[i] = divmod(N,b)
        return np.array(res).reshape(-1, 1)

