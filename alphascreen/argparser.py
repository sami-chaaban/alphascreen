import optparse
import sys
import alphascreen

def argparse():

    """
    optparse is used to initialize all command-line options.
    """
    
    parser = optparse.OptionParser(usage="Usage: %prog [options]",
        version=alphascreen.__version__)

    jobsetup_opts = optparse.OptionGroup(
        parser, 'Job setup')

    jobsetup_opts.add_option("--parse",
        action="store", dest="table", type="string", default="", metavar='file',
        help="Path to the excel file (.xlsx) or table (.txt).")

    jobsetup_opts.add_option("--fragment",
        action="store", dest="fraglen", type="string", default="", metavar='fragment-length',
        help="Approximate fragment length. Default is 500. For shorter fragments, try 250.")

    jobsetup_opts.add_option("--overlap",
        action="store", dest="overlap", type="int", default=50, metavar='overlap-length',
        help="Sequence is extended by this amount on either side of slices. Default is 50.")

    jobsetup_opts.add_option("--exhaustive",
        action="store_true", dest="exhaustive", default=False,
        help="Run every protein in column A against every protein in column B.")

    jobsetup_opts.add_option("--ignore_self",
        action="store_true", dest="ignoreself", default=False,
        help="Ignore any instances where both proteins are the same.")

    jobsetup_opts.add_option("--dimerize",
        action="store", dest="dimerize", type="string", default="", metavar='dimerize-id',
        help="Uniprot ID to dimerize. Alternatively, provide a text file (.txt) with a single column list of uniprot IDs to dimerize.")

    jobsetup_opts.add_option("--dimerize_all",
        action="store_true", dest="dimerize_all", default=False,
        help="Dimerize all proteins.")

    jobsetup_opts.add_option("--dimerize_all_except",
        action="store", dest="dimerize_except", type="string", default="", metavar='ids-not-to-dimerize',
        help="Dimerize all proteins except the uniprot ID provided here. Alternatively, provide a text file (.txt) with a single column list of uniprot IDs to NOT dimerize. Everything else will be dimerized.")

    jobsetup_opts.add_option("--consider",
        action="store", dest="consider", type="string", default="", metavar='sequence-to-consider',
        help="Uniprot ID and sequence range to consider. Example: \"Q86VS8/1/200\" only considers amino acids 1-200 for uniprot ID Q86VS8. Alternatively, provide a text file with a single column list of sequences to consider with a similar syntax (id/start/end).")

    jobsetup_opts.add_option("--alphafold_exec",
        action="store", dest="alphafold_exec", type="string", default="colabfold2", metavar='executable',
        help="Colabfold executable. Default is \"colabfold2\"")

    jobsetup_opts.add_option("--columnA",
        action="store", dest="columnA", type="string", default="SWISS-PROT Accessions Interactor A", metavar='columnA-name',
        help="Name of column heading for uniprot IDs for first interactors.")

    jobsetup_opts.add_option("--columnB",
        action="store", dest="columnB", type="string", default="SWISS-PROT Accessions Interactor B", metavar='columnB-name',
        help="Name of column heading for uniprot IDs for second interactors.")

    jobsetup_opts.add_option("--focus",
        action="store", dest="focus", type="string", default="", metavar='uniprot-id',
        help="Uniprot ID to focus on. This means that it will the first chain in any predictions that contain it.")

    jobsetup_opts.add_option("--dontwrite",
        action="store_true", dest="dontwrite", default=False,
        help="Pass if you don't want to write out files.")
    
    parser.add_option_group(jobsetup_opts)

    check_opts = optparse.OptionGroup(
        parser, 'Check runs')

    check_opts.add_option("--check",
        action="store_true", dest="check", default=False,
        help="Check how many jobs are left.")

    check_opts.add_option("--write_unfinished",
        action="store_true", dest="check_write", default=False,
        help="Check how many jobs are left and write a new script to only run those jobs.")

    parser.add_option_group(check_opts)

    analyze_opts = optparse.OptionGroup(
        parser, 'Analysis')

    analyze_opts.add_option("--show_top",
        action="store", dest="threshold", type="float", default=-1, metavar="threshold",
        help="Summarize results for top hits (above the threshold specified here) in tables and pdfs.")

    analyze_opts.add_option("--show_all",
        action="store_true", dest="showall", default=False,
        help="Summarize results for all hits in tables and pdfs.")

    analyze_opts.add_option("--write_table",
        action="store_true", dest="writetable", default=False,
        help="Write the full results as a table ranked by iptm.")

    jobsetup_opts.add_option("--rankby",
        action="store", dest="rankby", type="string", default="pae", metavar='ptm-iptm-pae',
        help="Score by which models are ranked (pae, ptm, or iptm). Default is pae. This is used for both choosing the best model in a prediction and then ranking the best models in the summaries.")
    
    analyze_opts.add_option("--overwrite",
        action="store_true", dest="overwrite", default=False,
        help="Overwrite snapshot pngs if they already exist. Default is False.")

    parser.add_option_group(analyze_opts)

    """
    The rest of the function parses the input and generates a dictionary for use in decisiontree.py
    """

    #Get the passed arguments from command-line
    options,args = parser.parse_args()

    #If there are less than 1 arguments passed, there were no options passed.
    if len(sys.argv) < 1:
            #parser.print_help()
            #print("\n>> You did not pass the arguments properly. For help, run \"starparser -h\".")
            print("\n>> Help: alphascreen -h\n")
            sys.exit()

    #Initialize an empty dictionary to place all the parameters in
    params={}

    #Place the passed parameters (or default values if none were passed) into the params dictionary
    for i in options.__dict__.items():
        params[i[0]] = i[1]
        
    #The dictionary is the main input to decisiontree.py
    return(params)