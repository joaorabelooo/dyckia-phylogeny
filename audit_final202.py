#!/usr/bin/env python3
"""Audit ML and PAUP SVD bootstrap trees using only the Python standard library."""

import csv
import re
import statistics
from datetime import datetime
from pathlib import Path

PROJECT = Path('/home/jrabelo/tese/filogenia_dyckia')
RUN = PROJECT / 'final202_trees.2jgge2'
OUTGROUP = 'ZD1_Deu_meziana'
EXCLUDED = {
    'Z42_D_magalhaesii', 'H39_D_reitzii_10',
    'D48_D_machrisiana_Leme_3291', 'H40_D_pseudodelicata',
    'D13_D_densiflora', 'D142_D_simulans',
}


def unquote(value):
    if value.startswith("'") and value.endswith("'"):
        return value[1:-1].replace("''", "'")
    return value


def parse_newick(text):
    tokens = re.findall(r"[(),:;]|'(?:[^']|'')*'|[^\s(),:;]+", text)
    position = 0

    def subtree():
        nonlocal position
        children = []
        if tokens[position] == '(':
            position += 1
            children.append(subtree())
            while tokens[position] == ',':
                position += 1
                children.append(subtree())
            if tokens[position] != ')':
                raise ValueError('Parêntese ausente no Newick')
            position += 1
        label = ''
        if tokens[position] not in {',', ')', ':', ';'}:
            label = unquote(tokens[position])
            position += 1
        length = None
        if tokens[position] == ':':
            length = float(tokens[position + 1])
            position += 2
        return {'children': children, 'label': label, 'length': length}

    root = subtree()
    if tokens[position:] != [';']:
        raise ValueError('Newick contém tokens adicionais ou falta ponto e vírgula')
    return root


def load_tree(path, nexus=False):
    text = path.read_text(encoding='utf-8')
    text = re.sub(r'\[.*?\]', '', text, flags=re.S)
    translation = {}
    if nexus:
        match = re.search(r'\btranslate\s+(.*?);', text, re.I | re.S)
        if match:
            for key, value in re.findall(
                r"(\d+)\s+('(?:[^']|'')*'|[^,\s]+)", match.group(1)
            ):
                translation[key] = unquote(value)
        trees = re.findall(r'\btree\s+[^=;]+=(.*?;)', text, re.I | re.S)
        if len(trees) != 1:
            raise ValueError(f'{path}: esperava uma árvore, encontrei {len(trees)}')
        text = trees[0]
    root = parse_newick(text.strip())
    leaves = []

    def visit(node):
        if node['children']:
            for child in node['children']:
                visit(child)
        else:
            node['label'] = translation.get(node['label'], node['label'])
            leaves.append(node['label'])

    visit(root)
    if len(leaves) != len(set(leaves)) or any(not label for label in leaves):
        raise ValueError(f'{path}: rótulos terminais vazios ou duplicados')
    return root, frozenset(leaves)


def collect_splits(root, universe, method):
    splits = {}

    def visit(node, is_root=False):
        if not node['children']:
            return frozenset([node['label']])
        side = frozenset().union(*(visit(child) for child in node['children']))
        if not is_root and min(len(side), len(universe - side)) >= 2:
            oriented = universe - side if OUTGROUP in side else side
            if method == 'ML':
                values = tuple(float(value) for value in node['label'].split('/'))
                if len(values) != 2:
                    raise ValueError(f'Esperava SH/UFBoot, encontrei {node["label"]}')
            else:
                values = (node['length'],)
            if any(value is None or not 0 <= value <= 100 for value in values):
                raise ValueError(f'Suporte inválido: {values}')
            if oriented in splits and splits[oriented] != values:
                raise ValueError('Suportes diferentes para a mesma bipartição')
            splits[oriented] = values
        return side

    if visit(root, True) != universe:
        raise ValueError('Terminais inconsistentes na árvore')
    return splits


def incompatible(a, b, universe):
    return all((a & b, a - b, b - a, universe - (a | b)))


def svd_newick(node, root=False):
    if node['children']:
        result = '(' + ','.join(svd_newick(c) for c in node['children']) + ')'
        if not root:
            if node['length'] is None:
                raise ValueError('Ramo interno sem bootstrap no SVD')
            result += f'{node["length"]:g}'
        return result
    return "'" + node['label'].replace("'", "''") + "'"


def write_table(out, name, fields, rows):
    with (out / name).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    samples = (PROJECT / '00_metadata/samples_final_202.txt').read_text(encoding='utf-8').split()
    universe = frozenset(samples)
    if len(samples) != 202 or len(universe) != 202:
        raise ValueError('Lista final não contém 202 IDs únicos')
    if OUTGROUP not in universe or EXCLUDED & universe:
        raise ValueError('Grupo externo ausente ou amostras excluídas na lista final')

    paths = {
        'ML': RUN / 'ml/dyckia_final202.mfp_asc.treefile',
        'SVD100': RUN / 'svdquartets/results/bootstrap100/dyckia_final202.mindist100bp.bootstrap100.tre',
        'SVD10000': RUN / 'svdquartets/results/bootstrap100/dyckia_final202.mindist10000bp.bootstrap100.tre',
    }
    roots, splits = {}, {}
    for method, path in paths.items():
        root, observed = load_tree(path, nexus=method != 'ML')
        if observed != universe:
            raise ValueError(
                f'{method}: ausentes={sorted(universe-observed)}, extras={sorted(observed-universe)}'
            )
        roots[method] = root
        splits[method] = collect_splits(root, universe, method)
        if len(splits[method]) > 199:
            raise ValueError(f'{method}: mais de 199 splits internos')
    if len(splits['ML']) != 199:
        raise ValueError('ML final com número inesperado de splits internos')

    # Use a new output directory to keep successive audits separate.
    out = RUN / ('audit_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir()
    lines = ['VALIDAÇÃO: mesmos 202 terminais únicos nas três árvores.',
             'Seis exclusões confirmadas. Z24 conserva seu ID histórico.',
             'ML: SH-aLRT/UFBoot. SVD: bootstrap padrão, escala 0–100.']
    strong_ml = {s for s, v in splits['ML'].items() if v[0] >= 80 and v[1] >= 95}
    lines.append(f'Ramos ML com SH>=80 e UF>=95: {len(strong_ml)}')
    summary = []
    for method in ('SVD100', 'SVD10000'):
        values = [v[0] for v in splits[method].values()]
        row = {'method': method, 'taxa': 202, 'internal_splits': len(values),
               'unresolved_of_199': 199-len(values),
               'mean_BS': round(statistics.mean(values), 3),
               'median_BS': statistics.median(values),
               'BS_lt50': sum(v < 50 for v in values)}
        for threshold in (50, 70, 95, 100):
            row[f'BS_ge{threshold}'] = sum(v >= threshold for v in values)
        summary.append(row)
        lines.append('\n' + method + ': ' + ', '.join(f'{k}={v}' for k, v in row.items()))
        (out / f'{method}.BS_labels.no_lengths.nwk').write_text(
            svd_newick(roots[method], root=True) + ';\n', encoding='utf-8'
        )
    write_table(out, 'support_summary.tsv', list(summary[0]), summary)

    all_sides = set().union(*(set(s) for s in splits.values()))
    ordered = sorted(all_sides, key=lambda s: (len(s), tuple(sorted(s))))
    split_rows = []
    stable_rows = []
    for side in ordered:
        ml = splits['ML'].get(side)
        b100 = splits['SVD100'].get(side)
        b10k = splits['SVD10000'].get(side)
        row = {'n': len(side), 'ML_SH': ml[0] if ml else '',
               'ML_UFBoot': ml[1] if ml else '',
               'SVD100_BS': b100[0] if b100 else '',
               'SVD10000_BS': b10k[0] if b10k else '',
               'members': ';'.join(sorted(side))}
        split_rows.append(row)
        if b100 and b10k and b100[0] >= 70 and b10k[0] >= 70:
            stable_rows.append(row)
    write_table(out, 'all_splits.tsv', list(split_rows[0]), split_rows)
    write_table(out, 'shared_SVD_BS70.tsv', list(split_rows[0]), stable_rows)
    lines.append(f'\nSplits com BS>=70 nos dois SVD: {len(stable_rows)}')
    lines.append('Amostras dos grupos estão em shared_SVD_BS70.tsv.')

    conflicts = []
    for left, right in (('ML', 'SVD100'), ('ML', 'SVD10000'), ('SVD100', 'SVD10000')):
        a, b = splits[left], splits[right]
        lines.append(f'\n{left} x {right}: RF={len(set(a)^set(b))}, compartilhados={len(set(a)&set(b))}')
        left_strong = strong_ml if left == 'ML' else {s for s, v in a.items() if v[0] >= 70}
        right_strong = {s for s, v in b.items() if v[0] >= 70}
        pairs = [(x, y) for x in left_strong for y in right_strong if incompatible(x, y, universe)]
        shared_strong = left_strong & right_strong
        lines.append(f'Compartilhados e sustentados nos dois: {len(shared_strong)}. Pares incompatíveis sustentados: {len(pairs)}')
        for x, y in sorted(pairs, key=lambda pair: (tuple(sorted(pair[0])), tuple(sorted(pair[1])))):
            conflicts.append({'method_A': left, 'method_B': right,
                              'support_A': '/'.join(f'{v:g}' for v in a[x]),
                              'support_B': '/'.join(f'{v:g}' for v in b[y]),
                              'members_A': ';'.join(sorted(x)),
                              'members_B': ';'.join(sorted(y))})
    write_table(out, 'supported_conflicts.tsv',
                ['method_A', 'method_B', 'support_A', 'support_B', 'members_A', 'members_B'], conflicts)
    lines.append(f'\nAUDITORIA CONCLUÍDA\nArquivos: {out}')
    report = '\n'.join(lines) + '\n'
    (out / 'comparison.log').write_text(report, encoding='utf-8')
    print(report)


if __name__ == '__main__':
    main()
