package chiralsoftware.moonphase;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.atan;
import static java.lang.Math.cos;
import static java.lang.Math.floor;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
import static java.lang.Math.toDegrees;
import static java.lang.Math.toRadians;
import java.time.LocalDate;
import java.time.LocalDateTime;
import static java.time.temporal.ChronoField.DAY_OF_MONTH;
import static java.time.temporal.ChronoField.MILLI_OF_DAY;
import static java.time.temporal.ChronoField.MONTH_OF_YEAR;
import static java.time.temporal.ChronoField.YEAR;
import static java.time.temporal.JulianFields.JULIAN_DAY;
import java.util.logging.Logger;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;

/**
 * Calculate moon phases.
 * For reference, here is a Julian date converter:
 * https://aa.usno.navy.mil/data/docs/JulianDate.php
 */
public final class MoonPhase {

    private static final Logger LOG = Logger.getLogger(MoonPhase.class.getName());
    
    private MoonPhase() { throw new RuntimeException(); } 

    // Astronomical constants.
    // 1980 January 0.0
    private static final double epoch = 2444238.5D;

    // Constants defining the Sun's apparent orbit.
    // ecliptic longitude of the Sun at epoch 1980.0
    private static final double elonge = 278.833540D;

    // ecliptic longitude of the Sun at perigee
    private static final double elongp = 282.596403D;

    // eccentricity of Earth's orbit
    private static final double eccent = 0.016718D;

    // semi-major axis of Earth's orbit, km
    private static final double sunsmax = 1.495985e8D;

    // sun's angular size, degrees, at semi-major axis distance
    private static final double sunangsiz = 0.533128D;

    // Elements of the Moon's orbit, epoch 1980.0.
    // moon's mean lonigitude at the epoch
    private static final double mmlong = 64.975464D;

    // mean longitude of the perigee at the epoch
    private static final double mmlongp = 349.383063D;

    // mean longitude of the node at the epoch
    // private static final double mlnode = 151.950429D;
    // inclination of the Moon's orbit
    // private static final double minc = 5.145396D;
    // eccentricity of the Moon's orbit
    private static final double mecc = 0.054900D;

    // moon's angular size at distance a from Earth
    private static final double mangsiz = 0.5181D;

    // semi-major axis of Moon's orbit in km
    private static final double msmax = 384401.0D;

    // parallax at distance a from Earth
    // private static final double mparallax = 0.9507D;
    // synodic month (new Moon to new Moon)
    private static final double synmonth = 29.53058868D;

    // base date for E. W. Brown's numbered series of lunations (1923 Jan 16)
    // private static final double lunatbase = 2423436.0D;
    // Properties of the Earth.
    // radius of Earth in kilometres
    // private static final double earthrad = 6378.16D;
    // Mathematical constants.
    private static final double EPSILON = 1E-6D;

    // Handy mathematical functions.
    // Fix angle.
    private static double fixangle(double a) {
        final double b = a - 360.0 * floor(a / 360.0D);
        // Can't use Math.IEEEremainder here because remainder differs
        // from modulus for negative numbers.
        return b;
    }

    // Sin from degrees.
    private static double dsin(double d) {
        return sin(toRadians(d));
    }

    // Cos from degrees.
    private static double dcos(double d) {
        return cos(toRadians(d));
    }

    /** Convert a Julian day into year / month / day . */
    private static void julianToYearMonthDay(double td, MutableInt yy, MutableInt mm, MutableInt dd) {
        final LocalDate ld = LocalDate.MIN.with(JULIAN_DAY, (long) floor(td + 0.5));
        yy.setValue(ld.get(YEAR));
        mm.setValue(ld.get(MONTH_OF_YEAR));
        dd.setValue(ld.get(DAY_OF_MONTH));
    }

    // meanphase - calculates mean phase of the Moon for a given base date
    // and desired phase:
    // 0.0 New Moon
    // 0.25 First quarter
    // 0.5 Full moon
    // 0.75 Last quarter
    // Beware!!! This routine returns meaningless results for any other
    // phase arguments. Don't attempt to generalise it without understanding
    // that the motion of the moon is far more complicated that this
    // calculation reveals.
    private static double meanphase(double sdate, double phase, MutableDouble usek) {
        final MutableInt yy = new MutableInt();
        final MutableInt mm = new MutableInt();
        final MutableInt dd = new MutableInt();
        double k, t, t2, t3, nt1;

        julianToYearMonthDay(sdate, yy, mm, dd);

        k = (yy.getValue() + ((mm.getValue() - 1) * (1.0 / 12.0)) - 1900) * 12.3685;

        // Time in Julian centuries from 1900 January 0.5.
        t = (sdate - 2415020.0) / 36525;
        t2 = t * t; // square for frequent use
        t3 = t2 * t; // cube for frequent use
        k = floor(k) + phase;
        usek.setValue(k);
        nt1 = 2415020.75933 + synmonth * k + 0.0001178 * t2 - 0.000000155 * t3
                + 0.00033 * dsin(166.56 + 132.87 * t - 0.009173 * t2);

        return nt1;
    }

    // truephase - given a K value used to determine the mean phase of the
    // new moon, and a phase selector (0.0, 0.25, 0.5, 0.75),
    // obtain the true, corrected phase time
    private static double truephase(double k, double phase) {
        final double t, t2, t3, m, mprime, f;
        double pt;
        boolean apcor = false;

        k += phase;
        /* add phase to new moon time */
        t = k / 1236.85;
        /*
						 * time in Julian centuries from 1900 January 0.5
         */
        t2 = t * t; /* square for frequent use */
        t3 = t2 * t; /* cube for frequent use */
        
        pt = 2415020.75933 /* mean time of phase */
                + synmonth * k + 0.0001178 * t2 - 0.000000155 * t3 + 0.00033
                * dsin(166.56 + 132.87 * t - 0.009173 * t2);

        m = 359.2242 /* Sun's mean anomaly */
                + 29.10535608 * k - 0.0000333 * t2 - 0.00000347 * t3;
        mprime = 306.0253 /* Moon's mean anomaly */
                + 385.81691806 * k + 0.0107306 * t2 + 0.00001236 * t3;
        f = 21.2964 /* Moon's argument of latitude */
                + 390.67050646 * k - 0.0016528 * t2 - 0.00000239 * t3;
        if ((phase < 0.01) || (abs(phase - 0.5) < 0.01)) {
            /* Corrections for New and Full Moon. */
            pt += (0.1734 - 0.000393 * t) * dsin(m) + 0.0021 * dsin(2 * m)
                    - 0.4068 * dsin(mprime) + 0.0161 * dsin(2 * mprime)
                    - 0.0004 * dsin(3 * mprime) + 0.0104 * dsin(2 * f) - 0.0051
                    * dsin(m + mprime) - 0.0074 * dsin(m - mprime) + 0.0004
                    * dsin(2 * f + m) - 0.0004 * dsin(2 * f - m) - 0.0006
                    * dsin(2 * f + mprime) + 0.0010 * dsin(2 * f - mprime)
                    + 0.0005 * dsin(m + 2 * mprime);
            apcor = true;
        } else if ((abs(phase - 0.25) < 0.01 || (abs(phase - 0.75) < 0.01))) {
            pt += (0.1721 - 0.0004 * t) * dsin(m) + 0.0021 * dsin(2 * m)
                    - 0.6280 * dsin(mprime) + 0.0089 * dsin(2 * mprime)
                    - 0.0004 * dsin(3 * mprime) + 0.0079 * dsin(2 * f) - 0.0119
                    * dsin(m + mprime) - 0.0047 * dsin(m - mprime) + 0.0003
                    * dsin(2 * f + m) - 0.0004 * dsin(2 * f - m) - 0.0006
                    * dsin(2 * f + mprime) + 0.0021 * dsin(2 * f - mprime)
                    + 0.0003 * dsin(m + 2 * mprime) + 0.0004
                    * dsin(m - 2 * mprime) - 0.0003 * dsin(2 * m + mprime);
            if (phase < 0.5)
                /* First quarter correction. */
                pt += 0.0028 - 0.0004 * dcos(m) + 0.0003 * dcos(mprime);
            else
                /* Last quarter correction. */
                pt += -0.0028 + 0.0004 * dcos(m) - 0.0003 * dcos(mprime);
            apcor = true;
        }
        if (!apcor)
            throw new InternalError("invalid phase selector");
        return pt;
    }

    // / Find time of phases of the moon which surround the current
    // date. Five phases are found, starting and ending with the
    // new moons which bound the current lunation.
    public static void phasehunt5(double sdate, double[] phases) {
        double adate, nt1, nt2;
        final MutableDouble k1 = new MutableDouble();
        final MutableDouble k2 = new MutableDouble();

        adate = sdate - 45;
        nt1 = meanphase(adate, 0.0, k1);
        for (;;) {
            adate += synmonth;
            nt2 = meanphase(adate, 0.0, k2);
            if (nt1 <= sdate && nt2 > sdate)
                break;
            nt1 = nt2;
            k1.setValue(k2);
        }
        phases[0] = truephase(k1.getValue(), 0.0);
        phases[1] = truephase(k1.getValue(), 0.25);
        phases[2] = truephase(k1.getValue(), 0.5);
        phases[3] = truephase(k1.getValue(), 0.75);
        phases[4] = truephase(k2.getValue(), 0.0);
    }

    // phasehunt2 - find time of phases of the moon which surround the current
    // date. Two phases are found.
    void phasehunt2(double sdate, double[] phases, double[] which) {
        double adate, nt1, nt2;
        final MutableDouble k1 = new MutableDouble();
        final MutableDouble k2 = new MutableDouble();

        adate = sdate - 45;
        nt1 = meanphase(adate, 0.0, k1);
        for (;;) {
            adate += synmonth;
            nt2 = meanphase(adate, 0.0, k2);
            if (nt1 <= sdate && nt2 > sdate)
                break;
            nt1 = nt2;
            k1.setValue(k2);
        }
        phases[0] = truephase(k1.getValue(), 0.0);
        which[0] = 0.0;
        phases[1] = truephase(k1.getValue(), 0.25);
        which[1] = 0.25;
        if (phases[1] <= sdate) {
            phases[0] = phases[1];
            which[0] = which[1];
            phases[1] = truephase(k1.getValue(), 0.5);
            which[1] = 0.5;
            if (phases[1] <= sdate) {
                phases[0] = phases[1];
                which[0] = which[1];
                phases[1] = truephase(k1.getValue(), 0.75);
                which[1] = 0.75;
                if (phases[1] <= sdate) {
                    phases[0] = phases[1];
                    which[0] = which[1];
                    phases[1] = truephase(k2.getValue(), 0.0);
                    which[1] = 0.0;
                }
            }
        }
    }

    // kepler - solve the equation of Kepler
    private static double kepler(double m, double ecc) {
        double e, delta;

        e = m = toRadians(m);
        do {
            delta = e - ecc * sin(e) - m;
            e -= delta / (1 - ecc * cos(e));
        } while (abs(delta) > EPSILON);
        return e;
    }

    /**
     * Calculate phase of moon as a fraction.
     *
     * @param pdate time for which the phase is requested, as from jtime()
     * @param pphaseR Ref for illuminated fraction of Moon's disk
     * @param moonAgeDays Ref for age of moon in days
     * @param distR Ref for distance in km from center of Earth
     * @param angdiaR Ref for angular diameter in degrees as seen from Earth
     * @param sudistR Ref for distance in km to Sun
     * @param suangdiaR Ref for Sun's angular diameter
     * @return terminator phase angle as a fraction of a full circle (i.e., 0 to
     * 1)
     */
    public static double phase(double pdate, MutableDouble pphaseR,
            MutableDouble moonAgeDays, MutableDouble distR, MutableDouble angdiaR,
            MutableDouble sudistR, MutableDouble suangdiaR) {
        final double day, N, M, lambdaSun, ml, MM, Ev, Ae, A3, MmP, mEc, A4, lP, V, lPP, moonAge, moonPhase, moonDist, moonDFrac, moonAng, F, sunDistance, sunAng;
        double Ec;

        // Calculation of the Sun's position.
        day = pdate - epoch; // date within epoch
        N = fixangle((360 / 365.2422) * day); // mean anomaly of the Sun
        M = fixangle(N + elonge - elongp); // convert from perigee co-ordinates
        // to epoch 1980.0
        Ec = kepler(M, eccent); // solve equation of Kepler
        Ec = sqrt((1 + eccent) / (1 - eccent)) * tan(Ec / 2);
        Ec = 2 * toDegrees(atan(Ec)); // true anomaly
        lambdaSun = fixangle(Ec + elongp); // Sun's geocentric ecliptic
        // longitude
        // Orbital distance factor.
        F = ((1 + eccent * cos(toRadians(Ec))) / (1 - eccent * eccent));
        sunDistance = sunsmax / F; // distance to Sun in km
        sunAng = F * sunangsiz; // Sun's angular size in degrees

        // Calculation of the Moon's position.
        // Moon's mean longitude.
        ml = fixangle(13.1763966 * day + mmlong);

        // Moon's mean anomaly.
        MM = fixangle(ml - 0.1114041 * day - mmlongp);

        // Evection.
        Ev = 1.2739 * Math.sin(toRadians(2 * (ml - lambdaSun) - MM));

        // Annual equation.
        Ae = 0.1858 * Math.sin(toRadians(M));

        // Correction term.
        A3 = 0.37 * Math.sin(toRadians(M));

        // Corrected anomaly.
        MmP = MM + Ev - Ae - A3;

        // Correction for the equation of the centre.
        mEc = 6.2886 * sin(toRadians(MmP));

        // Another correction term.
        A4 = 0.214 * sin(toRadians(2 * MmP));

        // Corrected longitude.
        lP = ml + Ev + mEc - Ae + A4;

        // Variation.
        V = 0.6583 * sin(toRadians(2 * (lP - lambdaSun)));

        // True longitude.
        lPP = lP + V;

        // Calculation of the phase of the Moon.
        // Age of the Moon in degrees.
        moonAge = lPP - lambdaSun;

        // Phase of the Moon.
        moonPhase = (1 - cos(toRadians(moonAge))) / 2;

        // Calculate distance of moon from the centre of the Earth.
        moonDist = (msmax * (1 - mecc * mecc))
                / (1 + mecc * cos(toRadians(MmP + mEc)));

        // Calculate Moon's angular diameter.
        moonDFrac = moonDist / msmax;
        moonAng = mangsiz / moonDFrac;

        pphaseR.setValue(moonPhase);
        moonAgeDays.setValue(synmonth * (fixangle(moonAge) / 360.0));
        distR.setValue(moonDist);
        angdiaR.setValue(moonAng);
        sudistR.setValue(sunDistance);
        suangdiaR.setValue(sunAng);
        return toRadians(fixangle(moonAge)) / (2*PI);
    }

    public static void main(String[] args) {
        LOG.info("Calculating some moon phases");
        final LocalDateTime local = LocalDateTime.now();
        final long julianLong = local.getLong(JULIAN_DAY);
        final double julianDouble = julianLong + local.getLong(MILLI_OF_DAY) / (1000f * 60 * 60 * 24);
        LOG.info("Julian long: " + julianLong + " and julian double: " + julianDouble);
        
        final MutableInt yy = new MutableInt();
        final MutableInt mm = new MutableInt();
        final MutableInt dd = new MutableInt();
        julianToYearMonthDay(julianDouble, yy, mm, dd) ;
        
        final MutableDouble pphaseR = new MutableDouble();

        final MutableDouble moonAge = new MutableDouble();
        final MutableDouble distR = new MutableDouble();
        final MutableDouble angdiaR = new MutableDouble();
        final MutableDouble sudistR = new MutableDouble();
        final MutableDouble suangdiaR = new MutableDouble();

        final double phase = phase(julianDouble, pphaseR,
                moonAge, distR, angdiaR,
                sudistR, suangdiaR);
        LOG.info("Ok it's done. Moon age in days = " + moonAge.getValue() + " and the phase is: " + phase);
        
        final double[] phases = new double[5];
        phasehunt5(julianDouble, phases);
        
        for(int i = 0; i < phases.length; i++) {
            julianToYearMonthDay(phases[i], yy, mm, dd);
            LOG.info("Phase: " + i + " at: " + yy + " " + mm + " " + dd + " (julian: " + phases[i] + ")");
        }
    }
}
