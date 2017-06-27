def cue(x, um, u, up, p, r):
    s = 1.0
    if um > up:
        s = -1.0
    y = s*x
    if y <= s*um: 
        cue = 0.0
    elif y >= s*up: 
        cue = y-s*u
    else:
        cue = p*(x-um)**2 + r*(x-um)**3
    return cue