"""
CODENAME:     PhyRe
DESCRIPTION:  
Copyright (c) 2009 Ronald R. Ferrucci, Federico Plazzi, and Marco Passamonti..
Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:
The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
"""


from optparse import OptionParser
from PhyRe import (
    initialize_parameters,
    read_population_file,
    read_sample_file,
    path_length,
    atd_mean,
    atd_variance,
    euler,
    print_results
)


def main():
    sample_file, pop_file, p, d1, d2, ci, batch, path_lengths, missing, out = initialize_parameters()

    population = {}
    population, taxon, coef, path_lengths_dict = read_population_file(pop_file)
    sample_data = read_sample_file(sample_file)

    if path_lengths == 'n':
        coef, pop_n, path_lengths = path_length(population)

    atd, taxon_n, taxon = atd_mean(sample_data, sample_data.keys())
    v_td = atd_variance(taxon_n, sample_data.keys(), atd)
    e_results = euler(sample_data, atd, taxon_n)

    results = {
        'atd': atd,
        'vtd': v_td,
        'euler': e_results,
        'N': taxon_n,
        'n': len(sample_data),
        'taxon': taxon
    }

    print_results(results, pop_n, path_lengths_dict)

if __name__ == "__main__":
    main()