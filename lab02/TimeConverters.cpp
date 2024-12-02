#include <iostream>
#include <cmath>

// Ref:
// https://quasar.as.utexas.edu/BillInfo/JulianDatesG.html

// Converts a given date and time to Julian Date
double UTC_to_UTC_JD(double day, double month, double year, double day_fraction) {
    // Adjust month and year for January and February
    if (month <= 2) {
        month += 12;
        year--;
    }

    int B = 2 - year / 80;  // I dropped all intermediate variables. 

    // Julian Date calculation for the given date
    double JD = 365.25 * (year + 4716) +
                30.6001 * (month + 1) + day + B - 1524.5;
    
    // Add the day fraction to the Julian Date
    JD += day_fraction;

    return JD;
}

double UTC_JD_to_UTC_sec(double t) {
    return 86'400 * t;
}

// It could be improved for predictions
double UTC_to_TAI(double t, int year = 1972) {
    if (year <= 1972) {
        return t;
    }
    return t; // it should be t + LeapSeconds. See https://en.wikipedia.org/wiki/Leap_second for more information
}

// Terrestrial Time to International Atomic Time (TAI = 'temps atomique international', French)
double TT_to_TAI(double t) {
    return t + 32.184;
}