"""Reproduce Figure S16 with the 2019 RSS code."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from urllib.request import urlopen

import numpy as np

root = Path(__file__).resolve().parent
SUSIER_COMMIT = '9837e1007808afb4380ecdc4270f57bd847824c3'
RSPARSEPRO_COMMIT = '5990896d58053a78d4260f011aa83949ae770c0a'


def verify_checksum(path, expected):
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Checksum mismatch: {path}")


def download(url, path, checksum):
    with urlopen(url, timeout=60) as response:
        path.write_bytes(response.read())
    verify_checksum(path, checksum)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--adsp', type=Path, default=root / 'FigS16.ADSP_R5_Chr22_LD.json')
parser.add_argument('--out', type=Path, default=root / 'output')
args = parser.parse_args()
out = args.out.resolve()
out.mkdir(parents=True, exist_ok=True)
verify_checksum(args.adsp, 'b7aabe109de9f7df6f4361e04647513a2cff89cf51241ceb05bb7914fea86445')
adsp = json.loads(args.adsp.read_text())
rng = np.random.default_rng(20261107)
cases = []
for kind in ['identity', 'adsp', 'low_rank']:
    for replicate in range(3):
        p = 12
        if kind == 'identity':
            R = np.eye(p)
        elif kind == 'adsp':
            R = np.asarray(adsp['blocks'][replicate]['R'])
        else:
            X = rng.normal(size=(6, p))
            X -= X.mean(axis=0)
            X /= np.sqrt((X * X).sum(axis=0))
            R = X.T @ X
        z = rng.normal(size=p) * 3
        for a in [0.1, 0.5, 1.0, 5.0]:
            cases.append({'case': f'{kind}_{replicate}_{a}', 'kind': kind,
                          'a': a, 'R': R.tolist(), 'z': z.tolist()})
with TemporaryDirectory() as directory:
    work = Path(directory)
    archive = work / 'susieR.tar.gz'
    download(f'https://codeload.github.com/stephenslab/susieR/tar.gz/{SUSIER_COMMIT}', archive,
             '6c887dc459c9ca4140560897c5ed06af7204e75c9755bb142e1bc9d343bee3f1')
    source = work / 'rsparsepro_ld.py'
    download(f'https://raw.githubusercontent.com/zhwm/RSparsePro_LD/{RSPARSEPRO_COMMIT}/src/rsparsepro_ld.py',
             source, '6480474383858f1033f788cac2e82adf288d4261456bf14df1fd5e30389e5ea8')
    spec = importlib.util.spec_from_file_location('upstream_rsparsepro', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    (work / 'inputs.json').write_text(json.dumps(cases))
    subprocess.run(
        ['Rscript', str(root / 'rss_fits.R'), str(root), str(work), SUSIER_COMMIT], check=True,
    )
    rss_results = {row['case']: row for row in json.loads((work / 'rss_fits.json').read_text())}

single_effect, variant_pip = [], []
maximum_alpha_difference = 0.0
for case in cases:
    R, z, a = np.asarray(case['R']), np.asarray(case['z']), case['a']
    reference = module.RSparsePro(len(z), 3, R, a)
    for _ in range(100):
        alpha_old = reference.gamma.copy()
        reference.infer_tilde_b(z)
        reference.infer_q_beta(R)
        if np.linalg.norm(reference.gamma - alpha_old) < 1e-5:
            break

    result = rss_results[case['case']]
    alpha_ref, pip_ref = reference.gamma, reference.get_PIP()
    rss_alpha = np.asarray(result['alpha']).reshape(3, len(z)).T
    internal_alpha = np.asarray(result['internal_alpha']).reshape(3, len(z)).T
    # Apply the same PIP summary as RSparsePro.
    pip_rss = rss_alpha.max(axis=1).round(4)
    pip_internal = internal_alpha.max(axis=1).round(4)
    for alpha, pip in [(rss_alpha, pip_rss), (internal_alpha, pip_internal)]:
        np.testing.assert_allclose(alpha, alpha_ref, rtol=0, atol=1e-10)
        maximum_alpha_difference = max(
            maximum_alpha_difference, float(np.max(np.abs(alpha - alpha_ref)))
        )
        np.testing.assert_array_equal(pip, pip_ref)

    for j in range(len(z)):
        variant = {'case': case['case'], 'kind': case['kind'], 'a': a, 'variant': j + 1}
        variant_pip.append({
            **variant, 'rss2019': float(pip_rss[j]), 'rsparsepro': float(pip_ref[j]),
        })
        for k in range(3):
            single_effect.append({
                **variant, 'effect': k + 1,
                'rss2019': float(rss_alpha[j, k]), 'rsparsepro': float(alpha_ref[j, k]),
            })

summary = {
    'cases': len(cases),
    'single_effect_probabilities': len(single_effect),
    'variant_PIPs': len(variant_pip),
    'maximum_alpha_difference': maximum_alpha_difference,
    'identical_PIP_cases': len(cases),
}
comparison = out / 'FigS16.RSparsePro_RSS_Lambda_Comparison.json'
comparison.write_text(json.dumps({
    'summary': summary, 'single_effect': single_effect, 'variant_pip': variant_pip,
}, indent=2) + '\n')
print(json.dumps(summary, indent=2))

subprocess.run(
    ['Rscript', str(root / 'plot_s16.R'), str(comparison), str(out)], check=True,
)
