class VCFWizard
  get '/' do
    template_render('main', {}, nil, :cache_type => :none)
  end
end
