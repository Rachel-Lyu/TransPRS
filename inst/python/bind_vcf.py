def bind_vcf(nameLis, out):
    fs = [open(fname, 'r') for fname in nameLis]
    fout = open(out, 'w')
    lines = []
    for i in range(len(nameLis)):
        while True:
            ts = next(fs[i])
            if not ts.startswith('##'):
                if i > 0:
                    lines.append(ts.split('\t', 9)[-1].strip())
                else:
                    lines.append(ts.strip())
                break
    print('\t'.join(lines), file = fout)

    for line0 in fs[0]:
        lines = [next(fs[i]).split('\t', 9)[-1].strip() if i > 0 else line0.strip() for i in range(len(nameLis))]
        print('\t'.join(lines), file = fout)
    for f in fs:
        f.close()
    fout.close()