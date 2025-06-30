pdbs = Dir.glob("*_fragments")
pdbs.each do |pdb|
       if(!File.exists?("#{pdb}/start.200.9mers") || !File.exists?("#{pdb}/start.200.3mers" )|| !File.exists?("#{pdb}/run.fold.boinc.job"))
        Dir.chdir(pdb)
        name = "hugh2020_#{File.basename(pdb)}"
        system("/home/brunette/scripts/prepare_fold_fasta.sh #{name}\n")
        #system("boinc_submit run.fold.boinc.job")
        Dir.chdir("..")
       end
end
