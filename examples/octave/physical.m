% LightPipes for Octave Optical Toolbox
% actually, this file is part of the olson-tools package
%** \file
%* A multitude of constants (mks where applicable) organized into a
%* bunch of namespaces.
%*
%* This file was ported from a C++ header that is freely available on usenet
%* servers.  The reference for this header is thus comp.lang.cpp/2004-10/2217
%* by E. Robert Tisdale (10/17/04). 
%*
%* There were a lot of other comments on usenet for improving this file,
%* including British spelling of various units and so on.  This version is the
%* unaltered version from E. R. Tisdale.
%*
%* */

%namespace physical {
%    namespace unit { % conversion factor
                     radian = 1.0;
                     radians = radian;
                     rad = radian;
                     steradian = 1.0;
                     steradians = steradian;
                     sr = steradian;
                     dollar = 1.0;
                     dollars = dollar;
                     cent = dollar/100.0;
                     cents = cent;
                     dozen = 12.0;
                     doz = dozen;
                     dz = dozen;
%        namespace bakers {
                  bakers.dozen = 13.0;
                  bakers.doz = bakers.dozen;
                  bakers.dz = bakers.dozen;
%        }
                     gross = 12.0*dozen;
                     gro = gross;
                     quire = 25.0;
                     quires = quire;
                     ream = 500.0;
                     reams = ream;
                     percent = 1.0/100.0;
                     proof = percent/2.0;
                     karat = 1.0/24.0;
                     karats = karat;
                     mole = 6.0221367e+23;
                     moles = mole;
                     mol = mole;
%                    pi = 3.14159265358979323846*radians;
%        namespace arc {
                     arc.degree = pi/180.0;
                     arc.degrees = arc.degree;
                     arc.minute = arc.degree/60.0;
                     arc.minutes = arc.minute;
                     arc.min = arc.minute;
                     arc.second = arc.minute/60.0;
                     arc.seconds = arc.second;
                     arc.sec = arc.second;
                     arc.grade = 0.9*arc.degrees;
                     arc.grades = arc.grade;
%            namespace centesimal {
            arc.centesimal.minute = arc.grade/100.0;
            arc.centesimal.minutes = arc.centesimal.minute;
            arc.centesimal.min = arc.centesimal.minute;
            arc.centesimal.second = arc.grade/10000.0;
            arc.centesimal.seconds = arc.centesimal.second;
            arc.centesimal.sec = arc.centesimal.second;
%              }
%            }

        % SI units (mks)
        % length
                     meter = 1.0;
                     meters = meter;
                     m = meter;
                     kilometer = 1000.0*meters;
                     kilometers = kilometer;
                     km = kilometer;
                     decimeter = meters/10.0;
                     decimeters = decimeter;
                     dm = decimeter;
                     centimeter = meters/100.0;
                     centimeters = centimeter;
                     cm = centimeter;
                     millimeter = meters/1000.0;
                     millimeters = millimeter;
                     mm = millimeter;
                     micron = meter/1000000.0;
                     microns = micron;
                     nanometer = meter/1000000000.0;
                     nanometers = nanometer;
                     nm = nanometer;
                     decinanometer = meter/10000000000.0;
                     decinanometers = decinanometer;
                     Angstrom = decinanometer;
                     Angstroms = Angstrom;
                     Xunit = 1.00202e-13*meters;
                     Xunits = Xunit;
                     Fermi = meter/1000000000000000.0;
                     Fermis = Fermi;
        % area
                     hectare = 10000.0*meter*meter;
                     hectares = hectare;
                     ha = hectare;
        % volume
                     stere = meter*meter*meter;
                     steres = stere;
                     liter = stere/1000.0;
                     liters = liter;
                     l = liter;
                     milliliter = stere/1000000.0;
                     milliliters = milliliter;
                     ml = milliliter;
                     cc = milliliter;
%        namespace displacement {
            displacement.ton = stere;
            displacement.tons = displacement.ton;
            displacement.t = displacement.ton;
%        }
        % mass
                     kilogram = 1.0;
                     kilograms = kilogram;
                     kg = kilogram;
                     quintal = 100.0*kilograms;
                     quintals = quintal;
                     doppelzentner = quintal;
                     doppelzentners = doppelzentner;
                     gram = kilograms/1000.0;
                     grams = gram;
                     g = gram;
                     milligram = kilogram/1000000.0;
                     milligrams = milligram;
                     mg = milligram;
%        namespace metric { % weight
                  metric.carat = gram/5.0;
                  metric.carats = metric.carat;
                  metric.ton = 1000.0*kilograms;
                  metric.tons = metric.ton;
                  metric.t = metric.ton;
%        }
        % time
                     second = 1.0;
                     seconds = second;
                     sec = second;
                     s = second;
                     millisecond = second/1000.0;
                     milliseconds = millisecond;
                     ms = millisecond;
                     microsecond = second/1000000.0;
                     microseconds = microsecond;
                     us = microsecond;
                     nanosecond = second/1000000000.0;
                     nanoseconds = nanosecond;
                     picosecond = second/1000000000000.0;
                     picoseconds = picosecond;
                     minute = 60.0*seconds;
                     minutes = minute;
%                    min = minute;
                     hour = 60.0*minutes;
                     hours = hour;
                     hr = hour;
                     day = 24.0*hours;
                     days = day;
                     da = day;
                     week = 7.0*days;
                     weeks = week;
                     fortnight = 2.0*weeks;
                     fortnights = fortnight;
                     year = 365.2421896698*days;
                     years = year;
                     yr = year;
                     month = year/12.0;
                     months = month;
                     mo = month;
                     decade = 10.0*years;
                     decades = decade;
                     century = 100.0*years;
                     centuries = century;
                     millenium = 1000.0*years;
                     millenia = millenium;
        % temporal frequency
                     Hertz = 1.0/second;
                     Hz = Hertz;
                     kiloHertz = 1000.0*Hertz;
                     kHz = kiloHertz;
                     megaHertz = 1000000*Hertz;
                     MHz = megaHertz;
                     gigaHertz = 1000000000.0*Hertz;
                     GHz = gigaHertz;
                     teraHertz = 1000000000000.0*Hertz;
                     THz = teraHertz;
        % spacial frequency
                     diopter = 1.0/meter;
                     diopters = diopter;
        % speed
                     kph = kilometers/hour;
        % radioactivity
                     Becquerel = 1.0/second;
                     Becquerels = Becquerel;
                     Bq = Becquerel;
                     Rutherford = 1000000.0*Becquerels;
                     Rutherfords = Rutherford;
                     Curie = 3.7e+10*Becquerels;
                     Curies = Curie;
                     Ci = Curie;
        % force
                     Newton = kilogram*meter/(second*second);
                     Newtons = Newton;
                     N = Newton;
                     dyne = Newton/100000.0;
                     dynes = dyne;
                     dyn = dyne;
        % pressure
                     Pascal = Newton/(meter*meter);
                     Pascals = Pascal;
                     Pa = Pascal;
                     Barie = Pascal/10.0;
                     Baries = Barie;
                     Barye = Barie;
                     Baryes = Barye;
                     pieze = 1000.0*Pascals;
                     piezes = pieze;
                     pz = pieze;
                     bar = 10000.0*Pascals;
                     bars = bar;
                     Torr = 133.3224*Pascals;
                     atmosphere = 760.0*Torr;
                     atmospheres = atmosphere;
                     atm = atmosphere;
        % energy
                     Joule = Newton*meter;
                     Joules = Joule;
                     J = Joule;
                     erg = Joule/10000000.0;
                     ergs = erg;
                     kiloWatthour = 3600000.0*Joules;
                     kiloWatthours = kiloWatthour;
                     kWh = kiloWatthour;
        % power
                     Watt = Joule/second;
                     Watts = Watt;
                     W = Watt;
                     kiloWatt = 1000.0*Watts;
                     kiloWatts = kiloWatt;
                     kW = kiloWatt;
                     megaWatt = 1000000.0*Watts;
                     megaWatts = megaWatt;
                     MW = megaWatt;
                     milliWatt = Watt/1000.0;
                     milliWatts = milliWatt;
                     mW = milliWatt;
                     microWatt = Watt/1000000.0;
                     microWatts = microWatt;
                     uW = microWatt;
%       namespace dose { % energy
                    dose.Gray = Joule/kilogram;
                    dose.Grays = dose.Gray;
                    dose.Gy = dose.Gray;
                    dose.Sievert = dose.Gray;
                    dose.Sieverts = dose.Sievert;
                    dose.rad = dose.Gray/100.0;
                    dose.rads = dose.rad;
                    dose.rd = dose.rad;
%       }
        % electrical current
                     Ampere = 1.0;
                     Amperes = Ampere;
                     A = Ampere;
                     Biot = 10.0*Amperes;
                     Biots = Biot;
        % electrical potential
                     Volt = Watt/Ampere;
                     Volts = Volt;
                     V = Volt;
        % electrical resistance
                     Ohm = Volt/Ampere;
                     Ohms = Ohm;
        % electrical conductance
                     mho = 1.0/Ohm;
                     mhos = mho;
                     Siemens = mho;
                     S = Siemens;
        % elsectrical charge
                     Coulomb = Ampere*second;
                     Coulombs = Coulomb;
                     C = Coulomb;
                     Franklin = 3.33564e-10*Coulombs;
                     Franklins = Franklin;
        % electrical capacity
                     Farad = Coulomb/Volt;
                     Farads = Farad;
                     F = Farad;
        % magnetic flux
                     Weber = Volt*second;
                     Webers = Weber;
                     Wb = Weber;
                     Maxwell = Weber/100000000.0;
                     Maxwells = Maxwell;
                     M = Maxwell;
        % magnetic field B
                     Tesla = Weber/(meter*meter);
                     Teslas = Tesla;
                     T = Tesla;
                     Gauss = Tesla/10000.0;
                     gamma = Tesla/1000000000.0;
        % magnetic field H
                     Oerstedt = 79.57747*Ampere/meter;
                     Oerstedts = Oerstedt;
                     Oe = Oerstedt;
        % magnetic inductivity
                     Henry = Weber/Ampere;
                     Henrys = Henry;
                     H = Henry;
                     milliHenry = Henry/1000.0;
                     milliHenrys = milliHenry;
                     mH = milliHenry;
        % temperature
                     Kelvin = 1.0;
                     Kelvins = Kelvin;
                     K = Kelvin;
                     milliKelvin = Kelvin / 1000.0;
                     milliKelvins = milliKelvin;
                     mK = milliKelvin;
                     microKelvin = Kelvin / 1000000.0;
                     microKelvins = microKelvin;
                     uK = microKelvin;
        % luminous intensity
                     candela = 1.0;
                     candelas = candela;
                     cd = candela;
                     apostilb = candelas/meter/meter;
                     apostilbs = apostilb;
                     nit = apostilb;
                     nits = nit;
                     skot = apostilb/1000.0;
                     skots = skot;
                     stilb = 10000.0*apostilbs;
                     stilbs = stilb;
                     Blondel = apostilb/pi;
                     Blondels = Blondel;
                     Lambert = 10000.0*Blondels;
                     Lamberts = Lambert;
        % light flux
                     lumen = candela*steradian;
                     lumens = lumen;
                     lm = lumen;
        % light intensity
                     lux = lumens/meter/meter;
                     luxes = lux;
                     luces = lux;
                     lx = lux;
                     nox = lux/1000.0;
                     phot = lumens/centimeter/centimeter;
                     phots = phot;
%       namespace equivalent {
%             equivalent.lux = unit.lux/pi;
              equivalent.lux = lux/pi;
              equivalent.luxes = equivalent.lux;
              equivalent.luces = equivalent.lux;
              equivalent.lx = equivalent.lux;
%             equivalent.lumen = unit.lumen/pi;
              equivalent.lumen = lumen/pi;
              equivalent.lumens = equivalent.lumen;
              equivalent.lm = equivalent.lumen;
              equivalent.phot = apostilb/pi;
              equivalent.phots = equivalent.phot;
%       }
        % acceleration
                     Galileo = centimeters/second/second;
                     Galileos = Galileo;
        % standard gavitational acceleration at sea level
                     gravity = 9.80665*meters/second/second;
        % mass
                     kilohyl = kilogram*gravity*second*second/meter;
                     kilohyls = kilohyl;
                     hyl = kilohyl/1000.0;
                     hyls = hyl;

        % English Units
        % length
                     inch = 0.0254*meters;
                     inches = inch;
                     in = inch;
                     mil = inch/1000.0;
                     mils = mil;
                     point = inch/72.27;
                     points = point;
                     pt = point;
                     bottommeasure = inch/40.0;
                     bottommeasures = bottommeasure;
                     line = inch/12.0;
                     lines = line;
                     pica = 12.0*points;
                     picas = pica;
                     barleycorn = inch/3.0;
                     barleycorns = barleycorn;
                     finger = 7.0*inches/8.0;
                     fingers = finger;
                     palm = 3.0*inches;
                     palms = palm;
                     hand = 4.0*inches;
                     hands = hand;
                     link = 7.92*inches;
                     links = link;
                     li = link;
                     span = 9.0*inches;
                     spans = span;
                     foot = 12.0*inches;
                     feet = foot;
                     ft = foot;
                     cubit = 18.0*inches;
                     cubits = cubit;
                     yard = 3.0*feet;
                     yards = yard;
                     yd = yard;
                     nail = yard/16.0;
                     nails = nail;
                     ell = 45.0*inches;
                     ells = ell;
                     pace = 5.0*feet;
                     paces = pace;
                     fathom = 6.0*feet;
                     fathoms = fathom;
                     fm = fathom;
                     rod = 198.0*inches;
                     rods = rod;
                     rd = rod;
                     pole = rod;
                     poles = pole;
                     p = pole;
                     perch = rod;
                     perches = perch;
                     rope = 20.0*feet;
                     ropes = rope;
                     bolt = 40.0*yards;
                     bolts = bolt;
                     chain = 4.0*rods;
                     chains = chain;
                     ch = chain;
%       namespace Gunters {
%                Gunters.chain = unit.chain;
                 Gunters.chain = chain;
                 Gunters.chains = Gunters.chain;
%       }
%       namespace engineers {
               engineers.link = foot;
               engineers.links = engineers.link;
               engineers.chain = 100.0*feet;
               engineers.chains = engineers.chain;
%       }
                     skein = 120*yards;
                     skeins = skein;
                     furlong = 220*yards;
                     furlongs = furlong;
                     spindle = 14400*yards;
                     spindles = spindle;
%       namespace US {
                      US.cable_length = 120.0*fathoms;
                      US.cable_lengths = US.cable_length;
%       }
%       namespace British {
                 British.cable_length = 100.0*fathoms;
                 British.cable_lengths = British.cable_length;
%       }
%       namespace statute {
                 statute.mile = 5280.0*feet;
                 statute.miles = statute.mile;
                 statute.mi = statute.mile;
                 statute.league = 3.0*statute.miles;
                 statute.leagues = statute.league;
%       }
%       namespace nautical {
                nautical.mile = 1852*meters;
                nautical.miles = nautical.mile;
                nautical.nm = nautical.mile;
                nautical.league = 3.0*nautical.miles;
                nautical.leagues = nautical.league;
%       }
%       namespace marine = nautical;
                  marine = nautical;
%       namespace geodetic {
                geodetic.foot = (1200.0/3937.0)*meters;
                geodetic.feet = geodetic.foot;
                geodetic.ft = geodetic.foot;
%       }
%       namespace geographical {
            geographical.mile = nautical.mile;
            geographical.miles = geographical.mile;
            geographical.mi = geographical.mile;
%       }
                     parasang = 3.5*statute.miles;
                     parasangs = parasang;
                     arpentcan = 27.52*statute.miles;
                     arpentcans = arpentcan;
                     arpentlin = 191.835*feet;
                     arpentlins = arpentlin;
                     astronomical_unit = 1.49597871e11*meters;
                     astronomical_units = astronomical_unit;
                     AU = astronomical_unit;
                     lightyear = 9.4605e15*meters;
                     lightyears = lightyear;
                     ly = lightyear;
                     parsec = AU*radians/arc.second;
                     parsecs = parsec;
                     pc = parsec;
        % area
                     barn = 1.0e-28*meter*meter;
                     barns = barn;
                     b = barn;
                     circular_inch = 0.25*pi*inch*inch;
                     circular_inches = circular_inch;
                     circular_mil = 0.25*pi*mil*mil;
                     circular_mils = circular_mil;
                     sabin = foot*foot;
                     sabins = sabin;
                     square = 100.0*sabin;
                     squares = square;
                     are = 100.0*meter*meter;
                     ares = are;
                     a = are;
                     rood = 40.0*rod*rod;
                     roods = rood;
                     ro = rood;
                     acre = 4.0*roods;
                     acres = acre;
                     section = statute.mile*statute.mile;
                     sections = section;
                     homestead = section/4.0;
                     homesteads = homestead;
                     township = 36.0*sections;
                     townships = township;
        % volume
                     minim = 6.161152e-8*(m*m*m);
                     minims = minim;
                     drop = 0.03*cc;
                     drops = drop;
                     teaspoon = 4.928922*cc;
                     teaspoons = teaspoon;
                     tablespoon = 3.0*teaspoons;
                     tablespoons = tablespoon;
%       namespace US {
%           namespace liquid {
                  US.liquid. dram = 60.0*minims;
                  US.liquid. drams = US.liquid.dram;
                  US.liquid. dr = US.liquid.dram;
                  US.liquid. ounce = 8.0*US.liquid.drams;
                  US.liquid. ounces = US.liquid.ounce;
                  US.liquid. oz = US.liquid.ounce;
                  US.liquid. gill = 4.0*US.liquid.ounces;
                  US.liquid. gills = US.liquid.gill;
                  US.liquid. gl = US.liquid.gill;
                  US.liquid. pint = 4.0*US.liquid.gills;
                  US.liquid. pints = US.liquid.pint;
                  US.liquid. pt = US.liquid.pint;
                  US.liquid. quart = 2.0*US.liquid.pints;
                  US.liquid. quarts = US.liquid.quart;
                  US.liquid. qt = US.liquid.quart;
                  US.liquid. magnum = 2.0*US.liquid.quarts;
                  US.liquid. magnums = US.liquid.magnum;
                  US.liquid. gallon = 4.0*US.liquid.quarts;
                  US.liquid. gallons = US.liquid.gallon;
                  US.liquid. gal = US.liquid.gallon;
%           }
%           namespace dry {
                  US.dry.    pint = 550.61047*cc;
                  US.dry.    pints = US.dry.pint;
                  US.dry.    pt = US.dry.pint;
                  US.dry.    quart = 2.0*US.dry.pints;
                  US.dry.    quarts = US.dry.quart;
                  US.dry.    qt = US.dry.quart;
%           }
                  US.    peck = 8.0*US.dry.quarts;
                  US.    pecks = US.peck;
                  US.    pk = US.peck;
                  US.    bushel = 4.0*US.pecks;
                  US.    bushels = US.bushel;
                  US.    bu = US.bushel;
                  US.    barrel = 31.5*US.liquid.gallons;
                  US.    barrels = US.barrel;
                  US.    bbl = US.barrel;
                  US.    bl = US.barrel;
%       }
%       namespace British {
%           namespace fluid {
            British.fluid. drachm = 60.0*minims;
            British.fluid. drachms = British.fluid.drachm;
            British.fluid. dr = British.fluid.drachm;
            British.fluid. ounce = 8.0*British.fluid.drachms;
            British.fluid. ounces = British.fluid.ounce;
            British.fluid. oz = British.fluid.ounce;
            British.fluid. gill = 5.0*British.fluid.ounces;
            British.fluid. gills = British.fluid.gill;
            British.fluid. gi = British.fluid.gill;
            British.fluid. pint = 4.0*British.fluid.gills;
            British.fluid. pints = British.fluid.pint;
            British.fluid. pt = British.fluid.pint;
            British.fluid. quart = 2.0*British.fluid.pints;
            British.fluid. quarts = British.fluid.quart;
            British.fluid. qt = British.fluid.quart;
            British.fluid. gallon = 4.0*British.fluid.quarts;
            British.fluid. gallons = British.fluid.gallon;
            British.fluid. gal = British.fluid.gallon;
%             }
            British.     peck = 2.0*British.fluid.gallons;
            British.     pecks = British.peck;
            British.     pk = British.peck;
            British.     bushel = 4.0*British.pecks;
            British.     bushels = British.bushel;
            British.     bu = British.bushel;
            British.     barrel = 36.0*British.fluid.gallons;
            British.     barrels = British.barrel;
            British.     bbl = British.barrel;
            British.     bl = British.barrel;
%       }
                     noggin = 2.0*US.liquid.ounces;
                     noggins = noggin;
                     cup = 8.0*US.liquid.ounces;
                     cups = cup;
                     fifth = US.liquid.gallon/5.0;
                     fifths = fifth;
                     jeroboam = 4.0*fifths;
                     jeroboams = jeroboam;
                     firkin = 9.0*US.liquid.gallons;
                     firkins = firkin;
                     kilderkin = 18.0*US.liquid.gallons;
                     kilderkins = kilderkin;
                     strike = 2.0*US.bushels;
                     strikes = strike;
                     sack = 3.0*US.bushels;
                     sacks = sack;
                     coomb = 4.0*US.bushels;
                     coombs = coomb;
                     seam = 8.0*US.bushels;
                     seams = seam;
                     wey = 40.0*US.bushels;
                     weys = wey;
                     last = 80.0*US.bushels;
                     lasts = last;
                     register_ton = 100.0*(ft*ft*ft);
                     register_tons = register_ton;
                     register_tn = register_ton;
                     cord = 128.0*(ft*ft*ft);
                     cords = cord;
                     cordfoot = cord;
                     cordfeet = cordfoot;
                     boardfoot = 144.0*inch*inch*inch;
                     boardfeet = boardfoot;
                     timberfoot = foot*foot*foot;
                     timberfeet = timberfoot;
                     hogshead = 2.0*US.barrels;
                     hogsheads = hogshead;
                     pipe = 4.0*US.barrels;
                     pipes = pipe;
                     tun = 8.0*US.barrels;
                     tuns = tun;
        % mass
                     grain = 0.06479891*grams;
                     grains = grain;
                     gr = grain;
                     pennyweight = 24.0*grains;
                     dwt = pennyweight;
%       namespace apothecary { % weight
              apothecary.scruple = 20.0*grains;
              apothecary.scruples = apothecary.scruple;
              apothecary.s = apothecary.scruple;
              apothecary.dram = 3.0*apothecary.scruples;
              apothecary.drams = apothecary.dram;
              apothecary.dr = apothecary.dram;
              apothecary.ounce = 8.0*apothecary.drams;
              apothecary.ounces = apothecary.ounce;
              apothecary.oz = apothecary.ounce;
              apothecary.pound = 12.0*apothecary.ounces;
              apothecary.pounds = apothecary.pound;
              apothecary.lb = apothecary.pound;
%       }
%       namespace troy = apothecary;
%       namespace ap = apothecary;
%       namespace t = apothecary;
                  troy = apothecary;
%                 ap = apothecary;
%                 t = apothecary;
%       namespace avoirdupois { % weight
             avoirdupois.pound = 7000.0*grains;
             avoirdupois.pounds = avoirdupois.pound;
             avoirdupois.lb = avoirdupois.pound;
             avoirdupois.ounce = avoirdupois.pound/16.0;
             avoirdupois.ounces = avoirdupois.ounce;
             avoirdupois.oz = avoirdupois.ounce;
             avoirdupois.dram = avoirdupois.ounce/16.0;
             avoirdupois.drams = avoirdupois.dram;
             avoirdupois.dr = avoirdupois.dram;
%       }
%       namespace avdp = avoirdupois;
%       namespace av = avoirdupois;
                  avdp = avoirdupois;
                  av = avoirdupois;
                     stone = 14.0*avoirdupois.pounds;
                     stones = stone;
                     st = stone;
%       namespace US { % short
                     US. hundredweight = 100.0*avoirdupois.pounds;
                     US. cwt = US.hundredweight;
                     US. quarter = US.hundredweight/4.0;
                     US. quarters = US.quarter;
                     US. qr = US.quarter;
                     US. ton = 20.0*US.hundredweight;
                     US. tons = US.ton;
                     US. tn = US.ton;
                     US. deadweight = US.ton;
%       }
%       namespace British { % long
                 British.hundredweight = 112.0*avoirdupois.pounds;
                 British.cwt = British.hundredweight;
                 British.quarter = British.hundredweight/4.0;
                 British.quarters = British.quarter;
                 British.qr = British.quarter;
                 British.ton = 20.0*British.hundredweight;
                 British.tons = British.ton;
                 British.tn = British.ton;
%       }
%       namespace English = British;
%       namespace Imperial = British;
                  English = British;
                  Imperial = British;
                     crith = 0.0906*grams;
                     criths = crith;
                     bag = 94.0*avoirdupois.pounds;
                     bags = bag;
                     cental = 100.0*avoirdupois.pounds;
                     centals = cental;
                     weymass = 252.0*avoirdupois.pounds;
        % rate
                     mgd = 1000000.0*US.liquid.gallons/day;
                     cfs = foot*foot*foot/second;
                     minersinch = 1.5*foot*foot*foot/minute;
                     mpg = statute.miles/US.liquid.gallon;
        % speed
                     mph = statute.miles/hour;
                     knot = nautical.miles/hour;
                     knots = knot;
%       namespace admiralty {
              admiralty. knot = 6980.0*feet/hour;
              admiralty. knots = admiralty.knot;
%       }
        % force
                     poundal = avdp.pound*foot/(second*second);
                     poundals = poundal;
                     pdl = poundal;
                     lbf = avoirdupois.pound*gravity;
        % pressure
                     psi = lbf/(inch*inch);
        % energy
                     calorie = 4.1868*Joules;
                     calories = calorie;
                     cal = calorie;
                     kilocalorie = 1000.0*calories;
                     kilocalories = kilocalorie;
                     kcal = kilocalorie;
                     Frigorie = kilocalorie;
                     Frigories = Frigorie;
                     Btu = 1055.06*Joules;
                     therm = 10000.0*Btu;
                     therms = therm;
                     thermie = 1000000.0*calories;
                     thermies = thermie;
        % power
                     horsepower = 735.49875*Watts;
                     HP = horsepower;
        % electrical current
                     Gilbert = 0.795775*Amperes;
                     Gilberts = Gilbert;
        % temperature
                     Rankin = 1.8*Kelvins;
                     Rankins = Rankin;
                     R = Rankin;
        % luminous intensity
                     candle = 1.02*candelas;
                     candles = candle;
%       namespace Hefner {
%                 Hefner.candle = 0.9*unit.candles;
                  Hefner.candle = 0.9*candles;
                  Hefner.candles = Hefner.candle;
%       }
        % light intensity
                     foot_candle = lumens/foot/foot;
                     foot_candles = foot_candle;
                     fc = foot_candle;
                     foot_Lambert = candelas/foot/foot/pi;
                     foot_Lamberts = foot_Lambert;
%       namespace equivalent {
%             equivalent.foot_candle = unit.foot_candle/pi;
              equivalent.foot_candle = foot_candle/pi;
              equivalent.foot_candles = equivalent.foot_candle;
              equivalent.fc = equivalent.foot_candle;
%       }
%   }
%   namespace units = unit;

%   namespace constant {
%       using namespace units;
        % speed of light
                     c = 2.99792458e8*meters/second;
        % speed of sound
                     Mach = 331.46*meters/second;
        % Planck constant
                     h = 6.6260755e-34*Joule*seconds;
                     h_bar = h/(2.0*pi);
        % standard gavitational acceleration at sea level
%                    g = units.gravity;
                     g = gravity;
        % electron charge
                     e = 1.60217733e-19*Coulombs;
        % elevtron Volt
                     eV = e*V;
                     keV = 1000.0*eV;
                     MeV = 1000000.0*eV;
                     GeV = 1000000000.0*eV;
                     Rydberg = 13.6054*eV;
                     Rydbergs = Rydberg;
        % electron mass
                     m_e = 9.1093897e-31*kilograms;
        % proton mass
                     m_P = 1.6726231e-27*kilograms;
        % deuteron mass
                     m_D = 1875.61339*MeV/(c*c);
        % unified atomic mass unit
                     atomic_mass_unit = 1.6605402e-27*kilograms;
                     atomic_mass_units = atomic_mass_unit;
                     amu = atomic_mass_unit;
                     Dalton = atomic_mass_unit;
                     Daltons = Dalton;
        % permittivity of free space
                     epsilon = 8.854187817e-12*Farads/meter;
        % permeability of free space
                     mu = 12.566370614e-7*Newtons/(A*A);
        % fine-structure constant
                     alpha = 1.0/137.0359895;
        % classical electron radius
                     r_e = 2.81794092e-15*meters;
        % electron Compton wavelength
                     lambda_bar = 3.86159323e-13*meters;
        % Bohr radius
                     a_0 = 0.529177249e-10*meters;
        % wavelength of 1 eV/c particle
                     lambda_1eV = 1.23984244e-6*meters;
        % Thomson cross section
                     sigma_0 = 0.66524616*barns;
        % Bohr magneton
                     mu_B = 5.78838263e-11*MeV/Tesla;
        % nuclear magneton
                     mu_N = 3.15245166e-14*MeV/Tesla;
        % electron cyclotron frequency/field
                     E_M_e = 1.75881962e11*C/kg*(rad/(s*T));
        % proton cyclotron frequency/field
                     E_M_P = 9.5788309e7*C/kg*(rad/(s*T));
        % gravitational constant
                     G = 6.67259e-11*m*m*m/(kg*s*s);
        % Avogadro's constant
                     N_A = 6.0221367e23;
        % Boltzmann constant
                     K_B = 1.380658e-23*Joules/Kelvin;
        % molar volume, ideal gas at standard temperature and pressure
                     V_molar = 2.897756e-3*meter*Kelvins;
        % Stefan-Boltzmann constant
                     sigma_SB = 5.67051e-8*W/(m*m*K*K*K*K);
        %
%   }
%   namespace constants = constant;
%} 
