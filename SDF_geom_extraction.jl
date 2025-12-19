using Rho2sdf
using Rho2sdf.DataImport

# Název úlohy pro výstupní soubory
taskName = "3D_2x1x1_4Legs"

# Import VTU souboru (automaticky detekuje typ elementů a hustoty)
vtu_file = "3D_2x1x1_4Legs-SIMP_results.vtu"
(X, IEN, rho) = import_vtu_mesh(vtu_file)

# Validace importovaných dat
if !validate_vtu_mesh(X, IEN, rho)
    error("Chyba při validaci importovaných dat!")
end

# Konfigurace nastavení
options = Rho2sdfOptions(
    # threshold_density = nothing,      # Automatická detekce prahu hustoty
    threshold_density = 0.5093,      # Value for isocotour (0, 1)
    sdf_grid_setup = :manual,        # Manuální zadání sítě (budete dotázáni)
    export_input_data = false,       # Neexportovat vstupní data
    export_nodal_densities = true,   # Exportovat nodální hustoty
    export_raw_sdf = true,          # Exportovat nevyhlazené SDF
    rbf_interp = true,              # Vyhlazování pomocí RBF interpolace
    rbf_grid = :same,               # Stejná síť pro RBF (lze změnit na :fine)
    remove_artifacts = true,         # Odstranit artefakty z SDF
    element_type = HEX8,             # Automaticky detekováno, ale lze zadat explicitně
)

# Spuštění hlavního procesu konverze
result = rho2sdf(taskName, X, IEN, rho, options = options)

# Rozbalení výsledků
(fine_sdf, fine_grid, sdf_grid, sdf_dists) = result
