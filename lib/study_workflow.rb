Workflow.require_workflow "KinMut"

module StudyWorkflow

  helper :appris_pis do
    @index ||= begin
                 url = "http://appris.bioinfo.cnio.es/download/data/appris_data.principal.homo_sapiens.tsv.gz"
                 tsv = TSV.open(url, :key_field => 1, :fields => [2], :type => :single)
                 tsv.fields = ["Ensembl Transcript ID"]
                 tsv.key_field = "Ensembl Gene ID"
                 tsv = tsv.swap_id "Ensembl Transcript ID", "Ensembl Protein ID", :identifiers => Organism.transcripts(organism), :persist => true
                 tsv.to_single
               end
  end

  helper :dbSNP_rs do
    @dbSNP_rs ||= DbSNP.mutations.tsv(:persist => true, :key_field => "Genomic Mutation", :fields => [0], :type => :single)
  end

  task :mutation_overview => :tsv do
    mutations = study.all_mutations
    tsv = TSV.setup(mutations, :key_field => "Genomic Mutation", :fields => [], :type => :double, :namespace => organism)

    mutation_samples = {}

    study.cohort.each do |genotype|
      name = genotype.jobname
      genotype.each do |mutation|
        mutation_samples[mutation] ||= []
        mutation_samples[mutation] << name
      end
    end

    tsv.add_field "Sample" do |mutation,values|
      mutation_samples[mutation]
    end

    log :genes, "Finding overlapping genes"
    mutation_genes = Misc.process_to_hash(mutations){|m| m.genes}
    tsv.add_field "Ensembl Gene ID" do |mutation,values|
      mutation_genes[mutation]
    end

    log :mutated_isoforms, "Inferring AA mutations in principal isoforms"
    mutation_principal_isoform = Misc.process_to_hash(mutations){|m| 
      all_mis = m.mutated_isoforms
      m_pis = m.genes.collect{|list| appris_pis.values_at *list}
      all_mis.zip(m_pis).collect{|mis,pis| 
        next if mis.nil? or pis.nil? or mis.empty? or pis.compact.empty?
        pi = pis.compact.first
        mis.select{|mi| mi.index(pi) == 0}.first
      }
    }

    log :damage, "Obtaining damage scores"
    all_mis = mutation_principal_isoform.values.flatten.compact.uniq

    isoform_scores = MutEval.job(:dbNSFP, study, :mutations => all_mis, :method => nil).run
    converted_scores = isoform_scores.fields.select{|f| f =~ /converted/ } 
    isoform_scores = isoform_scores.slice(converted_scores)
    isoform_MA_score = MutEval.job(:dbNSFP, study, :mutations => all_mis, :method => 'sift').run
    isoform_SIFT_score = MutEval.job(:dbNSFP, study, :mutations => all_mis, :method => 'mutation_assessor').run

    log :registering, "Registering results"
    tsv.add_field "Mutated Isoform" do |mutation,values|
      mi = mutation_principal_isoform[mutation]
      mi.nil? ? nil : [mi]
    end

    tsv.add_field "Change" do |mutation,values|
      mi = values["Mutated Isoform"].first
      mi ? [mi.split(":").last] : nil
    end

    tsv.add_field "Mutation Assessor" do |mutation,values|
      mi = values["Mutated Isoform"].first
      [isoform_MA_score[mi]]
    end

    tsv.add_field "SIFT" do |mutation,values|
      mi = values["Mutated Isoform"].first
      [isoform_SIFT_score[mi]]
    end

    tsv.add_field "Damage average" do |mutation,values|
      mi = values["Mutated Isoform"].first
      scores = isoform_scores[mi]

      avg = scores.nil? ? nil : Misc.mean(scores.flatten.compact.collect{|s| s.to_f}.reject{|s| s == -999})
      [avg]
    end
    tsv.add_field "SNP" do |mutation,values|
      [dbSNP_rs[mutation]]
    end

    tsv
  end

  dep :mutation_overview
  task :kinase_analysis => :tsv do
    tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Genomic Mutation", "Mutated isoform", "KinMut"])
    mutation_overview = step(:mutation_overview).load

    all_mis = mutation_overview.slice("Mutated Isoform").
      values.compact.flatten.uniq
    MutatedIsoform.setup(all_mis, organism)
    all_mis = all_mis.select_by(:consequence){|c| c == "MISS-SENSE"}

    mi2kinmut = Misc.process_to_hash(all_mis){|mis| 
      ddd mis
      kinmut = KinMut.job(:predict, study, :list => mis).run 
      fff kinmut
      exit
    }

    mutation_overview.through do |mutation, values|
      mis = values["Mutated Isoform"]
      mis
      tsv[mutation] = mi2kinmut
    end
  end
end
