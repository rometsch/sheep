import numpy as np

orbitalElementNames =  ["a", "e", "i", "AscendingNode",
                        "Pericenter", "TrueAnomaly"]

class Particle:
    pass

class Planet(Particle):

    def __init__(self, orbitalElements, auxInfo = {}):
        for elem in orbitalElementNames:
            try:
                setattr(self, elem, orbitalElements[elem])
            except KeyError:
                raise TypeError("Orbital Element '{}' not given!".format(elem))


        for key, data in auxInfo.items():
            setattr(self, key, data)

    def orbitalElements(self):
        # return a dict containing the orbital elements
        return { elem : getattr(self, elem) for elem in orbitalElementNames }

class System:

    def __init__(self):
        self.particles = []

    def add(self, p):
        if not isinstance(p, Particle):
            raise TypeError()
        self.particles.append(p)

class PlanetIni(System):
    # Parse planet config files as found in PlutoCuda and Fargo(3d)
    # they contain the orbital elements as whitespace separated rows

    def __init__(self, fpath, varnames = None):
        super().__init__()
        self.fpath = fpath
        self.varnames = varnames if varnames is not None else []
        self.parse()
        self.planets = self.particles

    def parse(self):
        self.header = ""
        with open(self.fpath, 'r') as f:
            for l in f:
                if l.strip()[0] != "#":
                    break
                self.header += l

        data = np.genfromtxt(self.fpath)
        if data.ndim == 1:
            data = [data]
        for n, pd in enumerate(data):
            orbitalElements = { v : e for v,e in zip(self.varnames, pd) \
                                if v in orbitalElementNames }
            auxInfo = { v : e for v,e in zip(self.varnames, pd) \
                        if v not in orbitalElementNames }
            self.add( Planet(orbitalElements, auxInfo = auxInfo) )

    def write(self, outpath):
        with open(outpath, 'w') as of:
            of.write(self.header)
            for p in self.planets:
                line = "   ".join(["{}".format(getattr(p,v)) for v in self.varnames])
                of.write(line + "\n")

class PlutoPlanetIni(PlanetIni):

    def __init__(self, fpath):
        varnames = ['m', 'a', 'e', 'i', 'AscendingNode'
                    , 'Pericenter', 'TrueAnomaly', 'Feedback']
        super().__init__(fpath, varnames = varnames)
