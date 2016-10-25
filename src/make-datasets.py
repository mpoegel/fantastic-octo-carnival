#! /bin/python
# Python v3.5.2

from bs4 import BeautifulSoup
from os.path import exists
import re


def parse_lat_long_cities(state):
    with open('data/raw/latlong-cities.htm', 'r') as fp:
        content = ''.join(fp.readlines())
    soup = BeautifulSoup(content, 'html.parser')
    pattern = state + '\n\n(.+\n)*'
    res = re.search(pattern, str(soup.pre))
    cities = res.group(0).split('\n')
    cities = cities[2:-1]
    output = []
    ids = []
    i = 0
    for city in cities:
        bits = city.split('  ')
        output.append('{0} {1} {2}'.format(i, bits[1].strip(), bits[2].strip()))
        name = bits[3].split(',')[0]
        ids.append((i, name))
        i += 1
    with open('data/processed/{0}-latlong.dat'.format(state), 'w+') as fp:
        fp.write('{0}\n'.format(len(output)))
        for city in output:
            fp.write(city + '\n')
    with open('data/processed/{0}-cities.dat'.format(state), 'w+') as fp:
        fp.write('{0}\n'.format(len(ids)))
        for i in range(len(ids)):
            fp.write('{0} {1}\n'.format(i, ids[i]))
    return ids


def parse_cities_pop(filename, state, cities):
    city_to_pop = dict()
    with open(filename, 'r') as fp:
        fp.readline()
        for line in fp:
            bits = line.split(',')
            city = bits[2].replace('"', '').replace('city', '').strip()
            pop = int(bits[-1].strip())
            city_to_pop[city] = pop
    pops = []
    for city in cities:
        if (city[1] in city_to_pop):
            pops.append((city[0], city_to_pop[city[1]]))
        else:
            print('WARNING - city not found: {0}'.format(city[1]))
    with open('data/processed/{0}-cities.dat'.format(state), 'w+') as fp:
        fp.write('{0}\n'.format(len(cities)))
        for city in cities:
            fp.write('{0} {1}\n'.format(city[0], city[1]))
    with open('data/processed/{0}-pops.dat'.format(state), 'w+') as fp:
        fp.write('{0} {1}\n'.format(len(pops), len(cities)))
        for pop in pops:
            fp.write('{0} {1}\n'.format(pop[0], pop[1]))
    return pops


def load_cities(filename, state):
    res = []
    lat_longs = []
    i = 0
    with open(filename, 'r') as fp:
        for line in fp:
            bits = line.split(',')
            lat_longs.append((i, bits[0].strip(), bits[1].strip()))
            city = bits[2].strip()
            res.append((i, city))
            i += 1
    with open('data/processed/{0}-latlong.dat'.format(state), 'w+') as fp:
        fp.write('{0}\n'.format(len(lat_longs)))
        for latlong in lat_longs:
            fp.write('{0} {1} {2}\n'.format(latlong[0], latlong[1], latlong[2]))
    return res


def main():
    if (exists('data/external/MN-cities.csv')):
        cities = load_cities('data/external/MN-cities.csv', 'Minnesota')
    else:
        cities = parse_lat_long_cities('Minnesota')
    parse_cities_pop('data/external/MN-pops.csv', 'Minnesota', cities)


if (__name__ == '__main__'):
    main()
