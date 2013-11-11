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
    tsv = TSV.setup({}, :key_field => "Genomic Mutation", :fields => [])
    header, line = header vcf

    fields = line.sub(/^#/,'').split(/\s+/)

    info_subfields = header["INFO"]
    if info_subfields
      fields.concat info_subfields.keys
      info_pos = fields.index "INFO"
    end

    tsv.fields = fields
    tsv.key_field = "Genomic Mutation"
    while line = vcf.gets
      chr, position, id, ref, alt, *rest = parts = line.split(/\s+/)
      if info_subfields
        info = rest[info_pos - 5]
        subfield_values = {}
        info.split(";").each{|p| k,v = p.split "="; subfield_values[k] = v || "TRUE"}
        ddd subfield_values
        extra_values = subfield_values.values_at *info_subfields.keys
        parts.concat extra_values
      end
      position, alt = Misc.correct_vcf_mutation(position.to_i, ref, alt)
      mutation = [chr, position.to_s, alt * ","] * ":"
      tsv[mutation] = parts
    end

    tsv
  end
end
  
