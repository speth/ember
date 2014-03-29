class TimeSeriesWriter(object):
    def __init__(self, solver, options):
        self.solver = solver
        self.options = options

    def __call__(self, name):
        print 'TimeSeriesWriter:', dir(self.solver)
