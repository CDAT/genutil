#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def ieq(a, b, r):
    return (abs((a-b)/r) < 1.e-5)