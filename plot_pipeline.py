import analysis as al

combine = True

analysis = al.AURKAnalysis(combine=combine, store_arrays=True)
analysis.plot_salt_bridges()

analysis.plot_hbonds_on_all_residues(use_reference=False)
analysis.plot_hbonds_on_all_residues(use_reference=True)

analysis.plot_hbonds_on_all_residues(use_reference=False, count_waters=True)
analysis.plot_hbonds_on_all_residues(use_reference=True, count_waters=True)

analysis.correlate_hbonds_bridges()

