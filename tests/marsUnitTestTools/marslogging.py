class marslogging():
    def __init__(self,fn):
        self.buf = []
        self.fn = fn
        return
    def setlog(self,fn):
        self.buf = []
        self.fn = fn
        return
    def info(self, s):
        self.buf = self.buf + [s + '\n']
        return
    def close(self):
        file = open(self.fn, 'w')
        file.writelines(self.buf)
        file.close()
        return