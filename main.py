import matplotlib.pyplot as plt
import datetime
import pandas as pd
import numpy as np
from astral.sun import sun
from astral.geocoder import database, lookup
import matplotlib.colors as colors
import matplotlib.cm as cmx
import astropy.time as time


def solve_kepler(orb_phase, e, n_iter=28):
    M = 2*np.pi*orb_phase
    E = M
    for _ in range(n_iter):
        E = E - (E - e*np.sin(E) - M)/(1-e*np.cos(E))
    l = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2)) + np.pi
    return l


def phase_on_date(jd, jd0, P):
    fi = (jd - jd0)/P - int((jd - jd0)/P)
    return fi


def distance(e, l):
    return (1 - e**2)/(1 + e*np.cos(l))


def plot_objects(ax, midnight_time, numb_of_days):
    t = time.Time(midnight_time, format='datetime')

    objects = pd.read_csv('objects.csv')
    numb_of_objects = len(objects)

    name_positions = []
    names = []
    for indx, obj in objects.iterrows():
        P, e, jd0 = obj['Period'], obj['e'], obj['JD0']

        phase = phase_on_date(t.jd, jd0, P)
        orb_ph = np.linspace(0, 1, 10000)
        l = solve_kepler(orb_ph, e)
        ang = (np.pi - abs(l - np.pi))/np.pi

        print(phase)
        ymin, ymax = 1 - (indx+1)/numb_of_objects, 1 - indx/numb_of_objects

        name_positions.append(.5*(ymax-ymin) + ymin)
        names.append(obj['Name'])

        ax.axhline(indx/numb_of_objects, c = 'k', linewidth=.5)

        numb_of_cycles = int(numb_of_days/P)
        one_cycle_xs = np.linspace(0, obj['Period'], 10000)
        if numb_of_cycles == 0:
            ax.fill_between(one_cycle_xs - phase*obj['Period'], ang*(ymax - ymin) + ymin,
                            y2=ymin, color='C3', alpha=0.9, linewidth=1)

            ax.fill_between(one_cycle_xs - phase*obj['Period']+obj['Period'], ang*(ymax - ymin) + ymin,
                            y2=ymin, color='C3', alpha=0.9, linewidth=1)
        else:
            for n_cycle in range(numb_of_cycles+2):
                ax.fill_between(one_cycle_xs - phase*obj['Period']+ n_cycle*obj['Period'],
                                ang*(ymax - ymin) + ymin, y2=ymin, color='C3', alpha=0.9, linewidth=1)

        ax.set_yticks(name_positions)
        ax.set_yticklabels(names)


def main():
    time_now = datetime.datetime.now()
    timezone_delta = datetime.timedelta(hours=-13)
    honolulu_utc_delta = datetime.timedelta(hours=-10)
    one_day_delta = datetime.timedelta(days=1)

    # Time now on Haleakala, Hawaii
    haleakala_time = time_now + timezone_delta
    haleakala_midnight_time = haleakala_time - datetime.timedelta(hours=haleakala_time.hour,
                                                                  minutes=haleakala_time.minute,
                                                                  seconds=haleakala_time.second,
                                                                  microseconds=haleakala_time.microsecond)

    fig, ax = plt.subplots(1, 1, figsize=(15, 5))

    numb_of_days = 14
    midnights = np.arange(0, numb_of_days+1, 1)

    ax.set_xlim(0, numb_of_days)
    ax.set_ylim(0, 1)

    ax.set_yticks([])
    ax.set_xticks(midnights)
    date_labels = [(haleakala_midnight_time + i*one_day_delta).strftime('%d.%m') for i in range(numb_of_days + 1)]
    ax.set_xticklabels(date_labels)

    plot_objects(ax, haleakala_midnight_time, numb_of_days)

    city = lookup("Honolulu", database())

    time_positions, time_labels = [], []

    for day_number in range(numb_of_days + 1):
        s = sun(city.observer, date=haleakala_midnight_time + day_number*one_day_delta, tzinfo=city.timezone)
        time_dawn = s['dawn']
        time_sunset = s['sunset']

        dawn_hour = int(time_dawn.strftime('%H'))
        dawn_minute = int(time_dawn.strftime('%M'))
        dawn_position = (dawn_hour*60 + dawn_minute)/1440
        time_positions.append(dawn_position + day_number)
        time_labels.append(time_dawn.strftime('%H:%M'))

        sunset_hour = int(time_sunset.strftime('%H'))
        sunset_minute = int(time_sunset.strftime('%M'))
        sunset_position = (sunset_hour*60 + sunset_minute)/1440
        time_positions.append(sunset_position + day_number)
        time_labels.append(time_sunset.strftime('%H:%M'))

        ax.axvspan(day_number + dawn_position, day_number + sunset_position, alpha=.5, hatch='//', edgecolor='k')
        ax.vlines(day_number + dawn_position, ymin=0, ymax=1, color='k', linewidth=.5, linestyle='dashed')
        ax.vlines(day_number + sunset_position, ymin=0, ymax=1, color='k', linewidth=.5, linestyle='dashed')

    ax2 = ax.twiny()
    ax2.set_xlim(0, numb_of_days)
    ax2.set_xticks(time_positions[:-2])
    ax2.set_xticklabels(time_labels[:-2])
    ax2.xaxis.set_ticks_position('bottom')

    ax.vlines(midnights, ymin=0, ymax=1, color='k', linewidth=.5, linestyle='dotted')
    ax.xaxis.set_ticks_position('top')
    plt.tight_layout()
    plt.savefig('ephemeris.pdf')
    plt.show()


if __name__ == '__main__':
    main()
