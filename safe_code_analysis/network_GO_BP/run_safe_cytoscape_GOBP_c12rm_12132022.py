# Author: Xiang Zhang
# Date created: 06/09/2021
# Date last modified: 12/13/2022

import sys
import os
from os.path import expanduser
from pathlib import Path

# Add path to folder containing the modified safepy library
sys.path.append(expanduser('~') + '/Documents/Research_CB/GIN/safepy_network_annotation/safepy/')

import safe

#%matplotlib inline

def safe_test_plot(input_network_file, input_attribute_file, output_dir, metric='shortpath_weighted_layout', neibor_radius=0.15, dist_thre=0.65, min_neighbor=10, enrichment_threshold=0.05, bool_mt=True):
    sf = safe.SAFE()
    
    sf.load_network(network_file=input_network_file, node_key_attribute="shared name")

    #sf.plot_network()

    sf.load_attributes(attribute_file=input_attribute_file)
    
    sf.define_neighborhoods(node_distance_metric=metric, neighborhood_radius=neibor_radius)
    
    sf.enrichment_threshold = enrichment_threshold # May test with different threshold

    sf.compute_pvalues(multiple_testing=bool_mt)

    # Test with different min size
    #sf.attribute_enrichment_min_size = min_neighbor

    sf.define_top_attributes()

    sf.define_domains(attribute_distance_threshold=dist_thre)

    sf.trim_domains()
    
    input_stem = Path(input_network_file).stem
    output_file_dir = os.path.join(output_dir, input_stem)
    if not os.path.exists(output_file_dir):
        os.makedirs(output_file_dir)

    # Prints output
    sf.plot_composite_network(show_each_domain=True, save_fig='{output_dir}/{input_file}_viz.pdf'.format(output_dir=output_dir, input_file=input_stem))
    sf.plot_composite_network(show_each_domain=False, save_fig='{output_dir}/top_{input_file}_viz.pdf'.format(output_dir=output_dir, input_file=input_stem))
    sf.print_output_files(output_dir=output_file_dir)
    sf.save_network(output_file=os.path.join(output_file_dir, input_stem + ".gpickle"))
    
    return sf

if __name__ == '__main__':
    # Expanded network
    input_file = expanduser('~') + 'expanded_c12rm_no507Var2_041_30_min010.cys'

    # Read in GO BP annotation attribute file
    input_attribute_file = expanduser('~') + 'go_p_matrix_0_500_with_name.txt'

    output_dir = expanduser('~') + '/no507Var2_expanded_'

    output_dir2 = output_dir + '30_041_GOBP_0_500_shortpath_weighted_layout_radius0_12_distThre0_65_enrichThre0_2_length_over_pcc'

    sf_bp = safe_test_plot(input_file, input_attribute_file, output_dir2, metric='shortpath_weighted_layout', neibor_radius=0.12, dist_thre=0.65, enrichment_threshold=0.2)
