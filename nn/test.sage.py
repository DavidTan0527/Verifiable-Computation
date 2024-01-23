

# This file was *autogenerated* from the file test.sage
from sage.all_cmdline import *   # import sage library

_sage_const_7 = Integer(7); _sage_const_0 = Integer(0); _sage_const_3 = Integer(3)
import random
import torch
from torch import Tensor

F = Zmod(_sage_const_7 )

# tensor1 = torch.randint(0, F.order(), (3,3))
arr1 = [F(random.randint(_sage_const_0 , F.order())) for _ in range(_sage_const_3 )]
tensor1 = Tensor(arr1)

arr2 = [F.random_element() for _ in range(_sage_const_3 )]
tensor2 = Tensor(arr2)

x = torch.matmul(tensor1, tensor2)
print(vector(arr1).dot_product(vector(arr2)))
print(tensor1, tensor2)
print(x)

