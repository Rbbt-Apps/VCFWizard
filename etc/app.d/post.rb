$LOAD_PATH.unshift 'lib'
require 'rbbt/workflow'
require 'rbbt/entity/study'
require 'vcf2table'

Workflow.require_workflow "Appris"
Workflow.require_workflow "MutEval"
Workflow.require_workflow "ICGC"

class VCFWizard
  class << self
    attr_accessor :vcf_files, :genotype_files, :studies
  end
  self.vcf_files = Rbbt.var.VCFWizard.vcf_files
  self.genotype_files = Rbbt.var.VCFWizard.genotype_files
  self.studies = Rbbt.var.VCFWizard.studies

  #{{{ VCF

  get '/new_vcf' do
    template_render('new_vcf', {}, nil, :cache_type => :none)
  end

  get '/remove_vcf/:name' do
    name = params[:name]
    remove_vcf(name)
    redirect to('/')
  end

  post '/new_vcf' do
    vcf = consume_parameter :vcf
    vcf_file = consume_parameter :vcf__param_file
    vcf = fix_input(:text, vcf, vcf_file)

    name = consume_parameter :name

    organism = consume_parameter :organism

    raise "No name" if name.nil? or name.empty?
    raise "No VCF" if vcf.nil? or vcf.empty?

    name = Misc.sanitize_filename name

    vcf = VCF2Table.open(StringIO.new(vcf))
    vcf.namespace = organism

    save_vcf(vcf, name)

    redirect to(File.join('/vcf/', name))
  end

  get '/vcf/:name' do
    name = params[:name]
    name = Misc.sanitize_filename name

    template_render('vcf', {:name => name}, 'VCF: ' << name, :cache_type => :asynchronous, :check => [VCFWizard.vcf_files[name]])
  end

  helpers do
    def remove_vcf(name)
      FileUtils.rm(VCFWizard.vcf_files[name])
    end

    def save_vcf(vcf, name)
      Open.write(VCFWizard.vcf_files[name], vcf.to_s)
    end

    def load_vcf(name)
      VCFWizard.vcf_files[name].tsv
    end
  end


  #{{{ Genotype

  post '/new_genotype/:name' do
    name = consume_parameter :name
    variants = consume_parameter :variants
    organism = consume_parameter :organism

    raise "No variants" if variants.nil? or variants.empty?

    save_genotype(name, variants, organism, true)
    
    redirect to('/')
  end

  get '/remove_genotype/:name' do
    name = params[:name]
    remove_genotype(name)
    redirect to('/')
  end

  get '/genotype/:name' do
    name = params[:name]
    name = Misc.sanitize_filename name

    template_render('genotype', {:name => name}, 'Genotype: ' << name, :cache_type => :asynchronous, :check => [VCFWizard.genotype_files[name]])
  end

  helpers do
    def remove_genotype(name)
      FileUtils.rm(VCFWizard.genotype_files[name])
    end

    def save_genotype(name, variants, organism, watson = true)
      variants = GenomicMutation.setup(variants.dup, name, organism, watson)
      Open.write(VCFWizard.genotype_files[name], Annotated.tsv(variants, :all).to_s)
    end

    def load_genotype(name)
      Annotated.load_tsv(VCFWizard.genotype_files[name].tsv)
    end

    def all_genotypes
      VCFWizard.genotype_files.glob('*').collect{|file| File.basename file}
    end
  end

  #{{{ Study

  get '/new_study' do
    template_render('new_study', {}, nil, :cache_type => :none)
  end

  post '/new_study' do
    name = consume_parameter :name
    genotypes = consume_parameter :genotypes

    genotypes = genotypes.select{|k,v| v.to_s == "true" }.keys
    raise "No genotypes specified" if genotypes.length == 0

    genotypes = genotypes.collect{|name| load_genotype name}
    organisms = genotypes.collect{|g| g.organism}.uniq

    raise "Not all genotypes are aligned to the same organism(/build version): #{organisms.inspect}" if organisms.length > 1

    create_study(name, genotypes, organisms.first)

    redirect to("/entity/Study/" << name)
  end

  get '/remove_study/:name' do
    name = params[:name]
    remove_study(name)
    redirect to('/')
  end
  helpers do
    def create_study(name, genotypes, organism)
      dir = Study.study_dir[name].find(:user).sub('{USER}', user || 'guest')
      Misc.in_dir(dir) do
        FileUtils.mkdir_p dir.genotypes
        genotypes.each do |genotype|
          name = genotype.jobname
          Open.write("genotypes/#{ name }", genotype * "\n")
        end
        info = {:watson => true, :organism => organism}
        Open.write('metadata.yaml', info.to_yaml)
      end
    end

    def remove_study(name)
      FileUtils.rm_rf Study.study_dir[name].find(:user).sub('{USER}', user || 'guest')
    end
  end

  get '/kinases' do
    template_render('kinases', {}, nil, :cache_type => :none)
  end
end

study = "Esophageal_Adenocarcinoma-UK"

path = Path.setup('', nil, nil, 
                  :global => ICGC.root.find["{PATH}"], 
                  :user => Rbbt.var.studies.find(:lib)["{USER}/{PATH}"], 
                  :local => Rbbt.var.studies.local.find(:lib)["{PATH}"])

Study.study_dir = path
