require 'rbbt-util'

module VCF2Table
  def self.header_lines(vcf)
    header_lines = []
    while line = vcf.gets
      if line =~ /^##/
        header_lines << line
      else
        return [header_lines, line]
      end
    end
    return [header_lines, line]
  end

  def self.header(vcf)
    lines, next_line = header_lines(vcf)

    header = {}
    lines.each do |line|
      if line =~ /^##([A-Z]+)=<ID=(.*),Number=(.*),Type=(.*),Description="(.*)">/
        field, id, number, type, description = $1, $2, $3, $4, $5
        subfield = {:numer => number, :type => type, :description => description}
        header[field] ||= {}
        header[field][id] = subfield
      else
      end
    end

    return [header, next_line]
  end

  def self.open(vcf)
    header, line = header vcf

    if line =~ /#/
      fields = line.sub(/^#/,'').split(/\s+/)
      line = vcf.gets
    else
      fields = %w(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample)[0..line.split(/\s+/).length-1]
    end

    tsv = TSV.setup({}, :key_field => "Genomic Mutation", :fields => fields)
    while line
      chr, position, id, ref, alt, *rest = parts = line.split(/\s+/)
      position, alt = Misc.correct_vcf_mutation(position.to_i, ref, alt)
      mutation = [chr, position.to_s, alt * ","] * ":"
      tsv[mutation] = parts
      line = vcf.gets
    end

    # Unfold FORMAT fields for each sample
    if format_pos = tsv.fields.index("FORMAT")
      format_fields = tsv[tsv.keys.first][format_pos].split(":")

      sample_fields = tsv.fields.values_at *(format_pos+1..tsv.fields.length-1).to_a

      format_fields.each_with_index do |ifield,i|
        sample_fields.each_with_index do |sfield,j|
          tsv.add_field "#{sfield}:#{ ifield }" do |mutation, values|
            values[format_pos+1+j].split(":")[i]
          end
        end
      end
    end

    # Unfold INFO fields
    if info_pos = tsv.fields.index("INFO")
      info_fields = tsv.values.collect{|v| v[info_pos].split(";").collect{|p| p.partition("=").first}}.compact.flatten.uniq.sort

      info_fields.each do |ifield|
        tsv.add_field "INFO:#{ ifield }" do |mutation, values|
          v = values[info_pos].split(";").select{|p| p.partition("=").first == ifield}.first
          if v
            field, sep, value = v.partition "="
            sep.empty? ? "true" : value
          else
            "false"
          end
        end
      end
    end

    tsv
  end
end
  
