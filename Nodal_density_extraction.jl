using Rho2sdf
using Rho2sdf.DataImport
using Dates

# Seznam všech úloh ke zpracování
taskNames = [
    "3D_2x1x1_4Legs_16tol_r2.0",
    "3D_2x1x1_4Legs_08tol_r2.0",
    "3D_2x1x1_4Legs_04tol_r2.0",
    "3D_2x1x1_4Legs_02tol_r2.0",
    "3D_2x1x1_4Legs_01tol_r2.0",
    "3D_2x1x1_MBB_16tol_r2.0",
    "3D_2x1x1_MBB_08tol_r2.0",
    "3D_2x1x1_MBB_04tol_r2.0",
    "3D_2x1x1_MBB_02tol_r2.0",
    "3D_2x1x1_MBB_01tol_r2.0",
]

# Funkce pro zápis logu
function write_log_file(taskName::String, elapsed_time::Float64)
    log_filename = "$(taskName)-Nodal_log.txt"
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    open(log_filename, "w") do io
        println(io, "==================================================")
        println(io, "SIMP TOPOLOGY OPTIMIZATION SUMMARY")
        println(io, "==================================================")
        println(io, "Task name:           $taskName")
        println(io, "Total time:          $(round(elapsed_time, digits=2)) s")
        println(io, "Generated:           $timestamp")
        println(io, "==================================================")
    end

    println("Log saved to: $log_filename")
end

# Zpracování všech úloh
for taskName in taskNames
    println("\n" * "="^60)
    println("Processing: $taskName")
    println("="^60)

    vtu_file = "$(taskName)-SIMP.vtu"

    # Kontrola existence souboru
    if !isfile(vtu_file)
        println("WARNING: File $vtu_file not found, skipping...")
        continue
    end

    # Měření času výpočtu
    start_time = time()

    # Import dat
    (X, IEN, rho) = import_vtu_mesh(vtu_file)

    # Validace
    if !validate_vtu_mesh(X, IEN, rho)
        println("ERROR: Validation failed for $taskName, skipping...")
        continue
    end

    # Výpočet nodal density a export
    (mesh, ρₙ) = compute_and_export_nodal_densities(taskName, X, IEN, rho)

    # Konec měření času
    elapsed_time = time() - start_time

    # Zápis logu
    write_log_file(taskName, elapsed_time)

    println("Completed: $taskName in $(round(elapsed_time, digits=2)) s")
end

println("\n" * "="^60)
println("All tasks completed!")
println("="^60)
