from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI
import logging

class CoverageDiscrepancyUtil:
    def __init__(self,config):
        self.config = config
        #instantiate class
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.gfu = GenomeFileUtil(self.callback_url)
        self.ru = ReadsUtils(self.callback_url)
        self.gaa = GenomeAnnotationAPI(self.callback_url)

    def run(self,params):
        # Get the read library as deinterleaved fastq files
 #       reads_filepath = self.get_reads_filepaths(params['reads_input'])
        genome_filepath = self.get_genome_filepath(params['genome_input_ref'])
        print(params)
        output = {}
        return output

    def get_reads_filepaths(self,input_ref):
        """
        This function needs a description
        :param input_ref:
        :return:
        """
        # Get the read library as deinterleaved fastq files
        reads_params = {'read_libraries': [input_ref],
                        'interleaved': 'false',
                        'gzipped': None
                        }
        reads = self.ru.download_reads(reads_params)['files']
        files = [reads[input_ref]['files']['fwd']]
        if reads[input_ref]['files']['rev']:
            files.append(reads[input_ref]['files']['rev'])
        print('running on files:')
        for f in files:
            print(f)
        return files


    # def _get_assembly_info(self, ref):
    #     ''' given a ref to an assembly or genome, figure out the assembly and return its info '''
    #     info = self.ws.get_object_info3({'objects': [{'ref': ref}]})['infos'][0]
    #     obj_type = info[2]
    #     if obj_type.startswith('KBaseGenomeAnnotations.Assembly') or obj_type.startswith('KBaseGenomes.ContigSet'):
    #         return {'info': info, 'ref': ref, 'genome_ref': None}
    #
    #     if obj_type.startswith('KBaseGenomes.Genome'):
    #         # we need to get the assembly for this genome
    #         ga = GenomeAnnotationAPI(self.service_wizard_url)
    #         assembly_ref = ga.get_assembly({'ref': ref})
    #         # using the path ensures we can access the assembly even if we don't have direct access
    #         ref_path = ref + ';' + assembly_ref
    #         info = self.ws.get_object_info3({'objects': [{'ref': ref_path}]})['infos'][0]
    #         return {'info': info, 'ref': ref_path, 'genome_ref': ref}
    #
    #     raise ValueError('Input object was not of type: Assembly, ContigSet or Genome.  Cannot get Bowtie2 Index.')


    # Shouldn't need this if get_assembly_info works
    # def determine_genome_or_assembly(self,genome_input_ref):
    #     pass

    def get_genome_filepath(self,genome_ref):
        genome = self.gaa.get_assembly({'ref': genome_ref})
        print(genome)


    def call_bowtie2(self,reads_filepath,):
        pass

        # print("Input parameters: " + pformat(params))
        # object_ref = params['object_ref']
        # object_info = self.ws_client.get_object_info_new({"objects": [{"ref": object_ref}],
        #                                                    "includeMetadata": 1})[0]
        # object_type = object_info[2]
        #
        # self.config['ctx'] = ctx
        # prokka_runner = ProkkaUtils(self.config)
        #
        # if "KBaseGenomeAnnotations.Assembly" in object_type:
        #     return [prokka_runner.annotate_assembly(params, object_info)]
        # elif "KBaseGenomes.Genome" in object_type:
        #     return [prokka_runner.annotate_genome(params)]
        # else:
        #     raise Exception("Unsupported type" + object_type)
        # #END annotate

    # Shouldn't need this:
    # def get_assembly_ref_from_genome(self, genome_ref, genome_object):
    #     """
    #     Given a Genome object, fetch the reference to its Assembly object on the workspace.
    #     Arguments:
    #       ref is a workspace reference ID in the form 'workspace_id/object_id/version'
    #       ws_obj download workspace object for the genome
    #     Returns a workspace reference to an assembly object
    #     """
    #     # Extract out the assembly reference from the workspace data
    #     ws_data = genome_object['data']
    #     assembly_ref = ws_data.get('contigset_ref') or ws_data.get('assembly_ref')
    #     if not assembly_ref:
    #         name = genome_object['info'][1]
    #         raise TypeError('The Genome ' + name + ' has no assembly or contigset references')
    #     # Return a reference path of `genome_ref;assembly_ref`
    #     return genome_ref + ';' + assembly_ref

