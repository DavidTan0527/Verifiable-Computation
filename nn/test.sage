import random
import torch
from torch import Tensor

F = Zmod(7)

# tensor1 = torch.randint(0, F.order(), (3,3))
arr1 = [F(random.randint(0, F.order())) for _ in range(3)]
tensor1 = Tensor(arr1)

arr2 = [F.random_element() for _ in range(3)]
tensor2 = Tensor(arr2)

x = torch.matmul(tensor1, tensor2)
print(vector(arr1).dot_product(vector(arr2)))
print(tensor1, tensor2)
print(x)

