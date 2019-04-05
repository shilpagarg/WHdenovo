echo ========First simulate Illumina reads=======
cd ../
python3 whdenovo.py simulate illumina test/mom1.fasta 2.5 test/illuminaSim

echo ========Then simulate PacBio reads=======
python3 whdenovo.py simulate pacbio \
        test/illuminaSim/child.het2.5.cov30_1.fq\
        30 \
        test/pacbioSim \
        test/illuminaSim/mom1.het2.5.cov30.fasta \
        test/illuminaSim/mom2.het2.5.cov30.fasta \
        test/illuminaSim/dad1.het2.5.cov30.fasta \
        test/illuminaSim/dad2.het2.5.cov30.fasta \
        test/illuminaSim/child1.het2.5.cov30.fasta \
        test/illuminaSim/child2.het2.5.cov30.fasta

echo =======Run assembly of one individual=======
python whdenovo.py assemble \
       --illumina1 test/illuminaSim/mom.het2.5.cov30_1.fq test/illuminaSim/dad.het2.5.cov30_1.fq test/illuminaSim/child.het2.5.cov30_1.fq \
       --illumina2 test/illuminaSim/mom.het2.5.cov30_2.fq test/illuminaSim/dad.het2.5.cov30_2.fq test/illuminaSim/child.het2.5.cov30_2.fq \
       --pacbio test/pacbioSim/pacbio_child.fasta

echo ======Run assembly of trioasm=======
python whdenovo.py assemble \
       --illumina1 test/illuminaSim/mom.het2.5.cov30_1.fq test/illuminaSim/dad.het2.5.cov30_1.fq test/illuminaSim/child.het2.5.cov30_1.fq \
       --illumina2 test/illuminaSim/mom.het2.5.cov30_2.fq test/illuminaSim/dad.het2.5.cov30_2.fq test/illuminaSim/child.het2.5.cov30_2.fq \
       --pacbio test/pacbioSim/pacbio_mom.fasta test/pacbioSim/pacbio_dad.fasta test/pacbioSim/pacbio_child.fasta \
       --ped ped.ped

echo =======validation=======
python whdenovo.py validate temp_${PID}_${TIME}/bc1 test/pacbioSim/pacbio_child.fasta