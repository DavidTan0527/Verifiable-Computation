
class Linear:
    def __init__(self, size_in, size_out, G = ZZ):
        self.size_in, self.size_out = size_in, size_out
        self.weights = random_matrix(G, size_out, size_in)
        self.bias = random_matrix(G, size_out, 1)

    def forward(self, x):
        if x.degree() != self.size_in:
            raise ValueError("Dimension does not match")
        return self.weights * x + self.bias


