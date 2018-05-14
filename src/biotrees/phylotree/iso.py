def equal(t1, t2):
    t1.sort()
    t2.sort()
    return t1 == t2


def iso(t1, t2):
    t1.shape().sort()
    t2.shape().sort()
    return t1.shape() == t2.shape()
