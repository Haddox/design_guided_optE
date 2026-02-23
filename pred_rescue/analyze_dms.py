import numpy as np
import pandas as pd

# --- Load data ---
lib1 = pd.read_csv(
    '/Users/hhaddox/Documents/baker_lab/design_guided_optE/rescue_dms/dms_data/dms_data_library_1.csv'
)
lib2 = pd.read_csv(
    '/Users/hhaddox/Documents/baker_lab/design_guided_optE/rescue_dms/dms_data/dms_data_library_2.csv'
)

# Rename lib2 stability column to match lib1
lib2 = lib2.rename(columns={'stabilityscore': 'stabilityscore_avg'})

# Filter to SSM rows only before concatenating
lib1_ssm = lib1[lib1['group'] == 'SSM'].copy()
lib2_ssm = lib2[lib2['group'] == 'SSM'].copy()

stability_scores_df = pd.concat([lib1_ssm, lib2_ssm], ignore_index=True)


# --- Compute delta_stabilityscore ---
def get_parent_stability(parent_name, parent_names, stabilityscore_metric):
    if parent_name not in parent_names:
        return np.nan
    else:
        return float(stability_scores_df[
            stability_scores_df['name'] == '{0}_wt'.format(parent_name)
        ][stabilityscore_metric].iloc[0])


ssm_data = stability_scores_df.copy()

parent_names = set(ssm_data['parent_name'])
ssm_data['parent_stabilityscore_avg'] = ssm_data.apply(
    lambda row: get_parent_stability(
        row['parent_name'], parent_names, 'stabilityscore_avg'
    ), axis=1
)
ssm_data['delta_stabilityscore'] = (
    ssm_data['stabilityscore_avg'] - ssm_data['parent_stabilityscore_avg']
)
ssm_data['delta_stabilityscore'] = round(ssm_data['delta_stabilityscore'], 1)


# --- Print results ---
rescue = ssm_data[ssm_data['delta_stabilityscore'] >= 1.0].copy()

all_designs = sorted(ssm_data['parent_name'].unique())
print(f"Designs with data ({len(all_designs)} total):")
for design in all_designs:
    n = sum(rescue['parent_name'] == design)
    print(f"  {design}: {n} rescue mutation(s)")


# --- Categorize rescue mutations by polarity ---
polar_aa = 'RHKDENQST'
rescue['wt_polarity'] = rescue['wt_aa'].apply(
    lambda x: 'polar' if x in polar_aa else 'nonpolar'
)
rescue['mut_polarity'] = rescue['mut_aa'].apply(
    lambda x: 'polar' if x in polar_aa else 'nonpolar'
)
polarity_dict = {
    key: []
    for key in ['wt_polarity', 'polar', 'nonpolar']
}
for wt_polarity in ['polar', 'nonpolar']:
    data_i = rescue[rescue['wt_polarity'] == wt_polarity]
    polarity_dict['wt_polarity'].append(wt_polarity)
    for mut_polarity in ['polar', 'nonpolar']:
        polarity_dict[mut_polarity].append(
            sum(data_i['mut_polarity'] == mut_polarity)
        )

polarity_df = pd.DataFrame(polarity_dict)
polarity_df.set_index('wt_polarity', inplace=True)

print("\n\nPolarity of rescue mutations (rows=wt, cols=mut):")
print(polarity_df)


# --- Classify nonpolar-to-nonpolar mutations by size change ---
np_data = rescue[
    (rescue['wt_polarity'] == 'nonpolar') &
    (rescue['mut_polarity'] == 'nonpolar')
].copy()


def classify_by_size(aa):
    if aa in ['A', 'G', 'P']:
        return 'small'
    elif aa in ['V', 'L', 'I', 'M']:
        return 'medium'
    elif aa in ['F', 'W', 'Y']:
        return 'large'
    else:
        raise ValueError(aa)


np_data['wt_size'] = np_data['wt_aa'].apply(classify_by_size)
np_data['mut_size'] = np_data['mut_aa'].apply(classify_by_size)


def classify_size_change(wt_size, mut_size):
    if wt_size == 'small':
        if mut_size == 'small':
            return 'size maintained'
        else:
            return 'small-to-large'
    elif wt_size == 'medium':
        if mut_size == 'small':
            return 'large-to-small'
        elif mut_size == 'medium':
            return 'size maintained'
        else:
            return 'small-to-large'
    else:
        assert wt_size == 'large'
        if mut_size == 'small':
            return 'large-to-small'
        elif mut_size == 'medium':
            return 'large-to-small'
        else:
            return 'size maintained'


np_data['size_change'] = np_data.apply(
    lambda row: classify_size_change(row['wt_size'], row['mut_size']),
    axis=1
)

print("\nSize change among nonpolar-to-nonpolar rescue mutations:")
print(np_data['size_change'].value_counts().to_string())


# --- Save rescue mutations dataframe ---
import os
os.makedirs('results', exist_ok=True)

rescue['size_change'] = np_data['size_change']  # NaN for non-nonpolar pairs

output_df = rescue[
    ['parent_name', 'wt_aa', 'site_n', 'mut_aa',
     'delta_stabilityscore', 'wt_polarity', 'mut_polarity', 'size_change']
].copy()
output_df['site_n'] = output_df['site_n'].astype(int)
output_df = output_df.sort_values(['parent_name', 'site_n']).reset_index(drop=True)

output_df.to_csv('results/rescue_mutations.csv', index=False)
print("\nSaved results/rescue_mutations.csv")
