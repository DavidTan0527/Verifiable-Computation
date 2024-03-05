from tqdm import tqdm

def bsgs_discrete_log(a, base):
    if a.curve() != base.curve():
        raise ValueError("Cannot do discrete log on elements of different curves")
    
    n = a.curve().order()
    m = ceil(sqrt(n))
    table = {}
    curr = 0
    for j in range(m):
    # for j in tqdm(range(m)):
        if curr not in table:
            table[curr] = j
        curr += base

    a_inv_m = -m * base
    assert(a_inv_m + m*base == 0)
    y = a
    for i in range(m):
        if y in table:
            return i * m + table[y]
        y = y + a_inv_m

    
