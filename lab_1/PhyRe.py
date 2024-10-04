import sys
from optparse import OptionParser
from re import match
from random import sample


sample_file = sys.argv[1]
del sys.argv[1]
pop_file = sys.argv[1]
del sys.argv[1]

p: int = 1000
d1: int = 10
d2: int = 70
ci: str = 'y'
batch: str = 'n'
path_lengths: str = 'n'
missing: str = 'n'

parser = OptionParser()

d1 = int(sys.argv[1])
del sys.argv[1]
d2 = int(sys.argv[1])
del sys.argv[1]

parser.add_option('-o')
parser.add_option('-p', type='int')
parser.add_option('-c')
parser.add_option('-b')
parser.add_option('-l')
parser.add_option('-m')

(options, args) = parser.parse_args()

missing = options.m if options.m else 'n'
out = options.o if options.o else sample_file.split('.')[0]
p = options.p if options.p else 1000
ci = options.c if options.c else 'y'
batch = options.b if options.b else 'n'
path_lengths = options.l if options.l else 'n'

sample = {}
population = {}
output: str = out + '.out'

o = open(output, 'a')

save_out = sys.stdout
sys.stdout = open(output, 'w')

if batch == 'y':
    files = []
else:
    files = [sample_file]

index = {}
taxon = {}
coef = {}
path_lengths_dict = {}

for line in open(sample_file):
    if batch == 'y':
        j = line.strip()
        files.append(j)
    else:
        break

duplicates = []

for line in open(pop_file):
    if match('Taxon:', line):
        x = line.split()
        x.remove('Taxon:')
        for i in x:
            taxon.append(i)
            j = x.index(i)
            index[i] = j + 1
        continue

    elif match('Coefficients:', line):
        x = line.split()
        x.remove('Coefficients:')
        x = map(eval, x)

        for t in taxon:
            i = taxon.index(t)
            coef[t] = sum(x[i:])
            path_lengths_dict[t] = x[i]

        continue

    line.strip()
    x = line.split()

    species = x[0]
    population[species] = {}

    if species in sample.keys():
        duplicates.append(species)
    else:
        sample[species] = {}
        population[species] = {}

    if missing == 'y':
        mtax = ''
        for t in taxon:
            if x[index[t]] == '/':
                sample[species][t] = mtax
            else:
                sample[species][t] = x[index[t]]
                mtax = x[index[t]]

            population[species][t] = sample[species][t]

    else:
        for t in taxon:
            sample[species][t] = x[index[t]]
            population[species][t] = sample[species][t]

if len(duplicates) > 0:
    print("Population master list contains duplicates:")
    for i in duplicates:
        print(i, '\n')

def path_length(population):
    taxon_n = {}
    x = {}
    for t in taxon:
        taxon[t] = {}
        x[t] = [population[i][t] for i in sample]

        if taxon.index(t) == 0:
            for i in set(x[t]):
                taxon[t][i] = x[t].count(i)
        else:
            for i in set(x[t]):
                if i not in x[taxon[taxon.index(t) - 1]]:
                    taxon[t][i] = x[t].count(i)

        taxon_n[t] = len(taxon[t])

    n = [float(len(taxon[t])) for t in taxon]
    n.insert(0, 1.0)
    raw = []
    for i in range(len(n) - 1):
        j = i + 1

        if n[i] > n[j]:
            c = 1
        else:
            c = (1 - n[i] / n[j])

        raw.append(c)

    s = sum(raw)
    adj_co = [i * 100 / s for i in raw]

    coef = {}
    path_lengths = {}
    for i in range(len(taxon)):
        t = taxon[i]
        coef[t] = sum(adj_co[i:])
        path_lengths[t] = adj_co[i]

    return coef, taxon_n, path_lengths

if path_lengths == 'n':
    coef, pop_n, path_lengths = path_length(population)
if path_lengths == 'y':
    xxx, pop_n, yyy = path_length(population)
    del xxx, yyy

def atd_mean(data: dict, sample: list) -> tuple:
    """Calculates the average taxonomic distinctness."""
    N = len(sample)
    taxon = {}
    taxon_n = {}
    av_td = 0
    n = 0

    for t in taxon:
        taxon[t] = {}
        x = [data[i][t] for i in sample]
        for i in set(x):
            taxon[t][i] = x.count(i)

    for t in taxon:
        taxon_n[t] = sum([taxon[t][i] * taxon[t][j] for i in taxon[t] for j in taxon[t] if i != j])
        n = taxon_n[t] - n
        av_td += n * coef[t]
        n = taxon_n[t]
    av_td /= (N * (N - 1))

    return av_td, taxon_n, taxon

def atd_variance(taxon_n: dict, sample: list, atd: float) -> float:
    """Calculates the variance of taxonomic distinctness."""
    v_td = 0
    N = 0
    n = 0

    for t in taxon:
        n = taxon_n[t] - n
        v_td += n * coef[t] ** 2
        n = taxon_n[t]

    N = len(sample)
    n = N * (N - 1)

    v_td = (v_td - ((atd * n) ** 2) / n) / n

    return v_td

def euler(data: dict, atd: float, taxon_n: dict) -> dict:
    """Calculates Euler index for the data."""
    sample = data.keys()
    n = len(sample)
    td_min = 0
    N = 0
    for t in taxon:
        k = len(taxon[t])
        td_min += coef[t] * (((k - 1) * (n - k + 1) * 2 + (k - 1) * (k - 2)) - N)
        N += ((k - 1) * (n - k + 1) * 2 + (k - 1) * (k - 2)) - N

    td_min /= (n * (n - 1))

    taxon.reverse()
    tax_max = {}
    taxon_n = {}
    
    for t in taxon:
        tax_max[t] = []
        if taxon.index(t) == 0:
            tax_max[t] = []
            for i in range(len(taxon[t])):
                tax_max[t].append([])
            for i in range(len(taxon[t])):
                tax_max[t][i] = [sample[j] for j in range(i, n, len(taxon[t]))]
        else:
            tax_max[t] = []
            for i in range(len(taxon[t])):
                tax_max[t].append([])
                s = taxon[taxon.index(t) - 1]
                tax = [tax_max[s][j] for j in range(i, len(taxon[s]), len(taxon[t]))]

                for j in tax:
                    tax_max[t][i] += j
        tax_max[t].reverse()

    taxon.reverse()
    td_max = 0
    n = 0
    N = len(sample)
    for t in taxon:
        taxon_n[t] = sum([len(tax_max[t][i]) * len(tax_max[t][j]) for i in range(len(tax_max[t])) for j in range(len(tax_max[t])) if i != j])
        n = taxon_n[t] - n
        td_max += n * coef[t]
        n = taxon_n[t]

    td_max /= (N * (N - 1))

    ei = (td_max - atd) / (td_max - td_min)

    e_results = {'EI': ei, 'TDmin': td_min, 'TDmax': td_max}
    return e_results

print("Output from Average Taxonomic Distinctness\n")

def sample(sample_file: str) -> dict:
    """Loads samples from a file. """
    sample = {}
    for line in open(sample_file):
        if match('Taxon:', line):
            continue
        elif match('Coefficients:', line):
            continue

        x = line.split()
        species = x[0]
        sample[species] = population[species]

    return sample

results = {}

for f in files:
    sample_data = sample(f)
    f = f.split('.')
    f = f[0]

    results[f] = {}

    samp = sample_data.keys()

    atd, taxon_n, taxon = atd_mean(sample_data, samp)
    v_td = atd_variance(taxon_n, samp, atd)
    e_results = euler(sample_data, atd, taxon_n)

    results[f]['atd'] = atd
    results[f]['vtd'] = v_td
    results[f]['euler'] = e_results
    results[f]['N'] = taxon_n
    results[f]['n'] = len(sample_data)
    results[f]['taxon'] = taxon

N = len(sample.keys())

def print_results() -> None:
    """Prints the analysis results to the screen."""
    print("Number of taxa and path lengths for each taxonomic level:")

    for t in taxon:
        print('%-10s\t%d\t%.4f' % (t, pop_n[t], path_lengths_dict[t]))
        n = 0

    print("\n")

    for f in results:
        print("---------------------------------------------------")
        print("Results for sample: ", f, '\n')
        print("Dimension for this sample is", results[f]['n'], '\n\n')
        print("Number of taxa and pairwise comparisons at each taxon level:")

        n = 0
        for t in taxon:
            N = results[f]['N'][t] - n
            print('%-10s\t%i\t%i' % (t, len(results[f]['taxon'][t]), N))
            n = results[f]['N'][t]

        print("""\nNumber of pairwise comparisons is for pairs that differ \
at each level excluding comparisons that differ at upper levels""")
        print("\n")

        print("Average taxonomic distinctness      = %.4f" % results[f]['atd'])
        print("Variation in taxonomic distinctness = %.4f" % results[f]['vtd'])
        print("Minimum taxonomic distinctness      = %.4f" % results[f]['euler']['TDmin'])
        print("Maximum taxonomic distinctness      = %.4f" % results[f]['euler']['TDmax'])
        print("von Euler's index of imbalance      = %.4f" % results[f]['euler']['EI'])
        print('\n')

print_results()
print("---------------------------------------------------")

sys.stdout = save_out
sys.stdout = sys.__stdout__

if ci == 'y':
    output = out.split('_')[0] + '_funnel.out'
    o = open(output, 'a')

    save_out = sys.stdout
    sys.stdout = open(output, 'w')
    print("""Confidence limits for average taxonomic distinctness and variation in taxonomic distinctness
limits are lower 95% limit for AvTD and upper 95% limit for VarTD
""")
    print("Number of permutations for confidence limits =", p, '\n')

    ci_array = []
    x = []
    c_array = []

    def funnel(p: int, d1: int, d2: int) -> None:
        """Calculates confidence intervals for taxonomic distinctness."""
        pop = population.keys()

        dims = []
        up = []
        lo = []
        means = []

        print("dimension AvTD05%   AvTDmean  AvTD95%   AvTDup    VarTDlow   VarTD05%   VarTDmean  VarTD95%")
        for d in range(d1, d2 + 1):
            x.append(d)
            av_td_ci = []
            var_td_ci = []
            for j in range(p):
                rsamp = sample(pop, d)

                atd, taxon_n, taxon = atd_mean(population, rsamp)
                av_td_ci.append(atd)
                v_td = atd_variance(taxon_n, rsamp, atd)
                var_td_ci.append(v_td)

            av_td_ci.sort()
            var_td_ci.sort()

            av_td = av_td_ci[int(.05 * p)], sum(av_td_ci) / p, av_td_ci[int(.95 * p)], max(av_td_ci)
            var_td = min(var_td_ci), var_td_ci[int(.05 * p)], sum(var_td_ci) / p, var_td_ci[int(.95 * p)]

            dims.append(d)
            ci_array.append(av_td[0])
            c_array.append(av_td[1])
            print('%i        %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f' %
                  (d, av_td[0], av_td[1], av_td[2], av_td[3], var_td[0], var_td[1], var_td[2], var_td[3]))

    funnel(p, d1, d2)

    sys.stdout = save_out
    sys.stdout = sys.__stdout__
