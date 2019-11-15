# rules:
.PHONY : run_all diff_all clean_all

default:
	
run_all: run_VP run_ROM run_TAL run_DNA run_NHP
	@echo "-------------------------------------------"
	@echo "New files generated."

diff_all: diff_VP diff_ROM diff_TAL diff_DNA diff_NHP
	@echo "-------------------------------------------"
	@echo "File comparison finished."

clean_all: clean_VP clean_TAL clean_ROM clean_DNA clean_NHP
	@echo "-------------------------------------------"
	@echo "Cleaning finished."

run_VP:
	../../src/aa_2_cg.jl --cgpdb --psf --log VP.pdb

diff_VP:
	diff VP_cg.gro VP_cg_ref.gro || exit 0
	diff VP_cg.pdb VP_cg_ref.pdb || exit 0
	diff VP_cg.psf VP_cg_ref.psf || exit 0
	diff VP_cg.itp VP_cg_ref.itp || exit 0
	diff VP_cg.log VP_cg_ref.log || exit 0

clean_VP:
	rm -f VP_cg.{log,itp,psf,pdb,gro}

run_ROM:
	../../src/aa_2_cg.jl --cgpdb --psf --log --force-field-protein Clementi ROM.pdb

diff_ROM:
	diff ROM_cg.gro ROM_cg_ref.gro || exit 0
	diff ROM_cg.pdb ROM_cg_ref.pdb || exit 0
	diff ROM_cg.psf ROM_cg_ref.psf || exit 0
	diff ROM_cg.itp ROM_cg_ref.itp || exit 0
	diff ROM_cg.log ROM_cg_ref.log || exit 0

clean_ROM:
	rm -f ROM_cg.{log,itp,psf,pdb,gro}

run_TAL:
	../../src/aa_2_cg.jl --cgpdb --psf --log TAL.pdb

diff_TAL:
	diff TAL_cg.gro TAL_cg_ref.gro || exit 0
	diff TAL_cg.pdb TAL_cg_ref.pdb || exit 0
	diff TAL_cg.psf TAL_cg_ref.psf || exit 0
	diff TAL_cg.itp TAL_cg_ref.itp || exit 0
	diff TAL_cg.log TAL_cg_ref.log || exit 0

clean_TAL:
	rm -f TAL_cg.{log,itp,psf,pdb,gro}

run_DNA:
	../../tools/DNA_3SPN.2C/bdna_3spn2c_gen.sh DNA.fasta

diff_DNA:
	diff DNA_x3dna_new_cg.gro DNA_cg_ref.gro || exit 0
	diff DNA_x3dna_new_cg.pdb DNA_cg_ref.pdb || exit 0
	diff DNA_x3dna_new_cg.psf DNA_cg_ref.psf || exit 0
	diff DNA_x3dna_new_cg.itp DNA_cg_ref.itp || exit 0

clean_DNA:
	rm -f DNA_x3dna*

run_TBP:
	../../src/aa_2_cg.jl --pwmcos --pwmcos-scale 2.0 --pwmcos-shift -3.14159 --pfm TBP.pfm TBP.pdb

diff_TBP:
	diff TBP_cg_pwmcos.itp_patch TBP_cg_ref.itp || exit 0

clean_TBP:
	rm -f TBP_cg_pwmcos.itp_patch

run_NHP:
	../../src/aa_2_cg.jl -f NHP_IDR.toml NHP.pdb

diff_NHP:
	diff NHP_cg.itp NHP_cg_ref.itp || exit 0

clean_NHP:
	rm -f NHP_cg.{itp,gro,top}
