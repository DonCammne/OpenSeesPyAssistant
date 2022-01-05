%  Module for units conversion to the default ones (m, N, kg, s)
%    Carmine Schipani, 2021

unit.m = 1;
unit.N = 1;
unit.kg = 1;
unit.s = 1;

% Distance
unit.mm = unit.m*1e-3;
unit.cm = unit.m*1e-2;
unit.dm = unit.m*1e-1;
unit.km = unit.m*1e3;
unit.inch = unit.m*0.0254;
unit.ft = unit.m*0.3048;
unit.mile = unit.m*1609.34;

% Area
unit.mm2 = unit.mm*unit.mm;
unit.cm2 = unit.cm*unit.cm;
unit.dm2 = unit.dm*unit.dm;
unit.m2 = unit.m*unit.m;
unit.inch2 = unit.inch*unit.inch;
unit.ft2 = unit.ft*unit.ft;

% Volume
unit.mm3 = unit.mm*unit.mm*unit.mm;
unit.cm3 = unit.cm*unit.cm*unit.cm;
unit.dm3 = unit.dm*unit.dm*unit.dm;
unit.m3 = unit.m*unit.m*unit.m;
unit.inch3 = unit.inch*unit.inch*unit.inch;
unit.ft3 = unit.ft*unit.ft*unit.ft;

% Moment of inertia
unit.mm4 = unit.mm3*unit.mm;
unit.cm4 = unit.cm3*unit.cm;
unit.dm4 = unit.dm3*unit.dm;
unit.m4 = unit.m3*unit.m;
unit.inch4 = unit.inch3*unit.inch;
unit.ft4 = unit.ft3*unit.ft;

% Force
unit.kN = unit.N*1e3;
unit.MN = unit.N*1e6;
unit.GN = unit.N*1e9;
unit.kip = unit.N*4448.2216;

% Moment
unit.Nm = unit.N*unit.m;
unit.kNm = unit.kN*unit.m;
unit.MNm = unit.MN*unit.m;
unit.Nmm = unit.N*unit.mm;
unit.kNmm = unit.kN*unit.mm;
unit.MNmm = unit.MN*unit.mm;

% Mass
unit.t = unit.kg*1e3;
unit.pound = unit.kg*0.45359237;

% Pressure/Stress
unit.Pa = unit.N/unit.m2;
unit.kPa = unit.Pa*1e3;
unit.MPa = unit.Pa*1e6;
unit.GPa = unit.Pa*1e9;
unit.psi = unit.Pa*6894.76;
unit.ksi = unit.psi*1000;

% Time
unit.min = unit.s*60;
unit.hours = unit.min*60;

%TODO: check with Andy

