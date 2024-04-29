from src.utils import cache_load_gene_summary

from src.app_interface import set_up_interface, file_selector
from src.app_configs_menu import process_configs_menu
from src.app_haplotype_plot import generate_haplotype_plot
from src.app_abacus_plot import generate_abacus_plot
from src.app_worldmap_plot import generate_worldmap_plot

def main():
    placeholder = set_up_interface()
    
    filename, gene_id_selected = file_selector(placeholder)
    
    df_haplotypes, df_join, background_ns_changes = cache_load_gene_summary(filename)

    min_samples, sample_count_mode = process_configs_menu(gene_id_selected, df_haplotypes, df_join)

    ns_changes, df_haplotypes_set = generate_haplotype_plot(df_haplotypes, gene_id_selected, background_ns_changes, min_samples, sample_count_mode)
    
    generate_abacus_plot(ns_changes, df_join, min_samples, df_haplotypes_set, gene_id_selected)
    
    generate_worldmap_plot(ns_changes, df_join, min_samples, gene_id_selected)

if __name__ == "__main__":
    main()