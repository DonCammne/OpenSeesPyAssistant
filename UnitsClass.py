# Module for units conversion to the default ones (m, N, kg, s)
#   Carmine Schipani, 2021

class Units():
    """Class with the units conversion and the definition of the units used as default (m, N, kg, s)
    """
    def __init__(self):
        # Units used in the library
        self.m = 1;
        self.N = 1;
        self.kg = 1;
        self.s = 1;

        # Distance
        self.mm = self.m*1e-3;
        self.cm = self.m*1e-2;
        self.dm = self.m*1e-1;
        self.km = self.m*1e3;
        self.inch = self.m*0.0254;
        self.ft = self.m*0.3048;
        self.mile = self.m*1609.34;

        # Area
        self.mm2 = self.mm*self.mm;
        self.cm2 = self.cm*self.cm;
        self.dm2 = self.dm*self.dm;
        self.m2 = self.m*self.m;
        self.inch2 = self.inch*self.inch;
        self.ft2 = self.ft*self.ft;

        # Volume
        self.mm3 = self.mm2*self.mm;
        self.cm3 = self.cm2*self.cm;
        self.dm3 = self.dm2**self.dm;
        self.m3 = self.m2*self.m;
        self.inch3 = self.inch2*self.inch;
        self.ft3 = self.ft2*self.ft;

        # Moment of inertia
        self.mm4 = self.mm3*self.mm;
        self.cm4 = self.cm3*self.cm;
        self.dm4 = self.dm3*self.dm;
        self.m4 = self.m3*self.m;
        self.inch4 = self.inch3*self.inch;
        self.ft4 = self.ft3*self.ft;

        # Force
        self.kN = self.N*1e3;
        self.MN = self.N*1e6;
        self.GN = self.N*1e9;
        self.kip = self.N*4448.2216;

        # Moment
        self.Nm = self.N*self.m;
        self.kNm = self.kN*self.m;
        self.MNm = self.MN*self.m;
        self.Nmm = self.N*self.mm;
        self.kNmm = self.kN*self.mm;
        self.MNmm = self.MN*self.mm;

        # Mass
        self.t = self.kg*1e3;
        self.pound = self.kg*0.45359237;

        # Pressure/Stress
        self.Pa = self.N/self.m2;
        self.kPa = self.Pa*1e3;
        self.MPa = self.Pa*1e6;
        self.GPa = self.Pa*1e9;
        self.psi = self.Pa*6894.76;
        self.ksi = self.psi*1000;

        # Time
        self.min = self.s*60;
        self.hours = self.min*60;



