#!/usr/bin/python

import subprocess
from collections import OrderedDict
import sys, os, math
import Bio.Phylo as bp
from Bio.Phylo import Newick
try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO



#-----------------------------------------------------------#
#           Classes used to parse centrifuge records        #
#-----------------------------------------------------------#


class CentrifugeRecord(object):
    """ This class used to represent centrifuge records in tabular result

        0      1       2     3         4           5        6
    readID uniqueID  taxID score secBestScore hitLength numMatches

    """
    __slots__=["readID", "uniqueID", "taxID", "score", "secBestScore",
               "hitLength", "numMatches"]

    record_dict = {} # {readID:CentrifugeRecord}

    def __init__(self, readID, uniqueID, taxID, score, secBestScore, hitLength, numMatches):
        self.readID = readID
        self.uniqueID = uniqueID
        self.taxID = taxID
        self.score = score
        self.secBestScore = secBestScore
        self.hitLength = hitLength
        self.numMatches = numMatches

        if not self.readID in CentrifugeRecord.record_dict:
            CentrifugeRecord.record_dict[self.readID] = [self]
        else:
            CentrifugeRecord.record_dict[self.readID].append(self)

    def __str__(self):
        return self.readID +"\t"+self.uniqueID+"\t"+\
               str(self.taxID)+"\t"+str(self.score)+"\t"+str(self.secBestScore)+\
               "\t"+str(self.hitLength)+"\t"+str(self.numMatches)




class CentrifugeRecordParser(object):
    """ This class used to parse centrifuge results iteratorly
    """
    def __init__(self, fh):
        """This parser is used to parse centrifuge result
        """
        self.handle = fh

    def _yield_records(self):
        """ This function used to yield mblastx record one by one
        """
        _file_handle = None
        # a opened file
        if hasattr(self.handle, 'read'):
            _file_handle = self.handle
        else:
            try:
                _file_handle = open(self.handle, "r")
            except Exception as e:
                print "cann't open file for read: %s !"%e

        for i, line in enumerate(_file_handle):
            if line.startswith("readID") or \
               line.startswith("\n") or \
               line.startswith(" ") or \
               line.startswith("\t"):
                continue
            line = line.strip().split("\t")
            #assert len(line) == 7, "The lenght of %d-th line is not 7 !!!"
            if len(line) != 7:
                print "Warning: the fields number of line %d is not 7, will be passed!!"%i
                print line
                continue
            yield CentrifugeRecord(*line)
    def __iter__(self):
        return self._yield_records()











#--------------------------------------------------------------#
#             BlastRecord and RecordParser                     #
#--------------------------------------------------------------#


class BlastRecord(object):
    """ This class used to represent tabular blast record in *.tab file, fields are

        0     1     2      3         4      5     6     7      8      9      10     11
    QID RefID Ident AlignLen MisMatch Gaps QStart QEnd RefStart RefEnd Evalue BitScore

    """

    __slots__=["QID", "RefID", "Ident", "AlignLen", "MisMatch", "Gaps",
               "QStart", "QEnd", "RefStart", "RefEnd", "Evalue", "BitScore"]



    def __init__(self,
                 QID=None, RefID=None,
                 Ident=None, AlignLen=None,
                 MisMatch=None, Gaps=None,
                 QStart=None, QEnd=None,
                 RefStart=None, RefEnd=None,
                 Evalue=None, BitScore=None
                 ):
        self.QID=QID
        self.RefID=RefID
        self.Ident=Ident
        self.AlignLen=AlignLen
        self.MisMatch=MisMatch
        self.Gaps=Gaps
        self.QStart=QStart
        self.QEnd=QEnd
        self.RefStart=RefStart
        self.RefEnd=RefEnd
        self.Evalue=Evalue
        self.BitScore=BitScore

    def __str__(self):
        """ string representation of blast record, the same with initial blast
            tabular output
        """
        _str=""
        _str += self.QID+"\t"
        _str += self.RefID+"\t"
        _str += self.Ident+"\t"
        _str += str(self.AlignLen)+"\t"
        _str += str(self.MisMatch)+"\t"
        _str += str(self.Gaps)+"\t"
        _str += str(self.QStart)+"\t"
        _str += str(self.QEnd)+"\t"
        _str += str(self.RefStart)+"\t"
        _str += str(self.RefEnd)+"\t"
        _str += str(self.Evalue)+"\t"
        _str += str(self.BitScore)
        return _str

class BlastRecordParser(object):
    """ This class used to parse tabular blast record in *.tab file, yield
        BlastRecord instance one by one
    """

    def __init__(self, fh):
        self.handle = fh


    def _yield_records(self):
        """ This function used to yield tabular blast records
        """
        _file_handle = None
        if not hasattr(self.handle, "read"):
            try:
                _file_handle = open(self.handle, "r")
            except Exception as e:
                print e
                print "can not open file for read !"
        else:
            _file_handle = self.handle

        for i, line in enumerate(_file_handle):
            if line.startswith("\n") or line.startswith("#"):
                continue
            line = line.strip().split("\t")
            assert len(line) == 12, "Bad blast record was found in %d-th line"%i+1
            yield BlastRecord(*line)

    def __iter__(self):
        return self._yield_records()





#-------------------------------------------------------------#
#    class of NCBI taxid, use to find lca, newick tree ...    #
#-------------------------------------------------------------#


class TaxIds(object):
    """ This class used to store taxid, sciName, parentTaxid and rank
        info, a path2file need as input to initialize. lines in path2file
        looks like below:

        taxid    scientificName    parentTaxid    rank\n
    """

    _taxdmp_ftp = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"

    def __init__(self, names=None, nodes=None):
        self.names = names
        self.nodes = nodes
        self.container = self._update_taxid_dict()


    def _update_taxid_dict(self, dir=os.getcwd()):
        """
        :param dir: tmp directory to put downloaded taxdmp.zip, and unzipped folder
        :return: taxid_dict
        """
        _taxid_dict = {} # {taxid: [scientificName, parent, rank]}

        # open names and nodes file for read, if given
        if self.names and self.nodes:
            try:
                names = open(self.names, "r")
                nodes = open(self.nodes, "r")
            except Exception as e:
                print "Cannot open names or nodes file for read: %s"%e

        # download taxdmp file, unzip, then open names and nodes file for read
        else:

            if os.path.exists(os.path.join(dir, "taxdmp.zip")):
                print "\n\nWARNING: taxdmp.zip exists, will be removed!\n\n"
                os.system("rm %s"%os.path.join(dir, "taxdmp.zip"))

            # download
            try:
                proc = subprocess.Popen(["wget",
                                               self._taxdmp_ftp, # ftp address for taxdmp
                                               "-P", dir    # directory to save taxdmp
                                             ])
                proc.wait()
            except Exception as e:
                print "Cannot download taxdmp file because: %s"%e
            # unzip
            try:
                proc = subprocess.Popen(["unzip",
                                         "-o",                              # overwrite files without prompting
                                         "-d", os.path.join(dir, "taxdmp"),  # extract files into directory
                                         os.path.join(dir, "taxdmp.zip")])
                proc.wait()
            except Exception as e:
                print "Error raised in unzipping taxdmp.zip file: %s"%e
            # read
            try:
                names = open(os.path.join(dir, "taxdmp", "names.dmp"), "r")
                nodes = open(os.path.join(dir, "taxdmp", "nodes.dmp"), "r")
            except Exception as e:
                print "Cannot open names or nodes file for read: %s"%e

        # read names.dmp,
        for line in names:
            if line.startswith("\n"):
                continue
            line = line.rstrip("\t|\n").split("\t|\t")
            if line[3] == "scientific name":
                # first get taxid:[scitificname]
                _taxid_dict.update({int(line[0]):[line[1]]})

        for line in nodes:
            if line.startswith("\n"):
                continue
            line = line.rstrip("\t|\n").split("\t|\t")
            taxid = int(line[0])
            parent = int(line[1])
            rank = line[2]
            if taxid in _taxid_dict:
                _taxid_dict[taxid].append(parent)
                _taxid_dict[taxid].append(rank)
            else:
                _taxid_dict.update({taxid:["None", parent, rank]})

        if not names.closed:
            names.close()
        if not nodes.closed:
            nodes.close()

        return _taxid_dict


    def get_parent(self, taxid):
        """ given a taxid, return its parent's taxid
        :param taxid: input query taxid
        :return:      the parent taxid of input query taxid
        """
        if not taxid:
            return None

        if not type(taxid) == int:
            try:
                taxid = int(str(taxid).strip())
            except Exception as e:
                print "taxid %s cannot converted to integer: %s !!"%(taxid, e)

        if self.container.has_key(taxid):
            return self.container[taxid][1]

    def get_sciName(self, taxid):
        """ given a taxid, return its scientific name
        :param taxid: input query taxid
        :return:      scientific name of input query taxid
        """
        if taxid is None:
            return None

        if not type(taxid) == int:
            try:
                taxid = int(str(taxid).strip())
            except Exception as e:
                print "taxid %s cannot converted to integer: %s !!"%(taxid, e)

        if self.container.has_key(taxid):
            return self.container[taxid][0]

    def get_rank(self, taxid):
        """ given a taxid, return its rank info
        :param taxid:  input query taxid
        :return:       the taxonomy rank of input query taxid
        """
        if not taxid:
            return None

        if not type(taxid) == int:
            try:
                taxid = int(str(taxid).strip())
            except Exception as e:
                print "taxid %s cannot converted to integer: %s !!"%(taxid, e)

        if self.container.has_key(taxid):
            return self.container[taxid][2]


    def get_path(self, taxid, toStr=False):
        """ given a taxid, return its path from root to this node
        :param taxid:  input query taxid
        :param toStr:  if True, return path in ";" separated string format, else in list format
        :return:       the taxonomic path from root to current taxid
        """
        if not taxid:
            return None

        if not type(taxid) == int:
            try:
                taxid = int(str(taxid).strip())
            except Exception as e:
                print "taxid %s cannot converted to integer: %s !!"%(taxid, e)

        _path = []
        _path.append(taxid)
        while taxid:
            taxid = self.get_parent(taxid)
            if taxid:
                _path.append(taxid)

            if taxid == 1:
                break
        # reverse list, making it from root to leaves
        _path.reverse()

        if toStr:
            _str = [self.get_sciName(taxid) for taxid in _path]
            return ";".join(_str)+";"
        else:
            return _path

    def get_query_rank_from_path(self, Taxid_Path, query_rank, toStr=False):
        """
        :param Taxid_Path: the path of taxid from root to leaves
        :param query_rank: the rank level want to query
        :param toStr:      if true, return scientific name, else return taxid
        :return:           the queried rank in scientific name or taxid
        """

        # assume only one meet the query_rank in one Taxid_Path
        ret = False
        for taxid in Taxid_Path:
            if not type(taxid) == int:
                try:
                    taxid = int(str(taxid).strip())
                except Exception as e:
                    print "taxid %s cannot converted to integer: %s !!"%(taxid, e)
            _rank = self.get_rank(taxid)
            if _rank == query_rank:
                ret = True
                if toStr:
                    return self.get_sciName(taxid)
                else:
                    return taxid
        # if no query_rank found in taxid_path, then return None
        if not ret:
            return None


    def _get_first_common_ancestor(self, taxid1, taxid2=None):
        """ this function take two taxid nodes as input, find their first
            common ancestor. This function should be used inside self.get_lca(),
            and should not be called outside.
        :param taxid1:  input first query taxid
        :param taxid2:  input second query taxid
        :return:        the common ancestor of these two taxid
        """
        # if two taxids are given, compare their pathes to find common_ancestor
        if taxid2:
            _path1 = self.get_path(taxid1)
            _path2 = self.get_path(taxid2)

            # find shorter path and longer path
            shorter, longer = (_path1, _path2) if len(_path1) <= len(_path2) \
                                else (_path2, _path1)

            if len(shorter) == 0:
                return longer.pop()

            # iterate over shorter path to check whether item in longer path
            common_ancestor = shorter.pop()
            while common_ancestor not in longer:
                try:
                    common_ancestor = shorter.pop()
                except IndexError:
                    print "%s or %s is not tracing back to the root !!!"%(
                            str(_path1), str(_path2))
                    break
            return common_ancestor

        # if only one taxid is given, then the renturn the last node in its path
        else:
            _path1 = self.get_path(taxid1)
            return _path1.pop()


    def get_lca(self, ListOfTaxid):
        """ given a list of taxid, return their lowest common ancestor
        :param ListOfTaxid:
        :return:
        """
        if len(ListOfTaxid) > 1:
            ListOfTaxid = [int(item) for item in ListOfTaxid]
            return reduce(self._get_first_common_ancestor, ListOfTaxid)

        elif len(ListOfTaxid) == 1:
            return self._get_first_common_ancestor(ListOfTaxid[0])

        else:
            print "\n\n No Taxid was given in ListOfTaxid !!\n\n"
            sys.exit(0)

    def path2newick(self, path2pathFile, node_fmt="taxid", out_fmt="newick"):
        """ This function take taxonomic path file as input, path should be consist
            of taxonomic id, not scitific name, because some scientific name are the
            same in different rank, but ids are unique.

            node_fmt = taxid / sciName

            out_fmt = newick / phyloxml ...

        """
        path, fileName = os.path.split(path2pathFile)
        basename = os.path.splitext(fileName)[0]
        outFile = os.path.join(path, basename+"2tree_"+node_fmt+"."+out_fmt)


        with open(path2pathFile, "r") as pathFile:

            # read in pathFile, and store node info into nodes
            nodes = {} # data format {"node_name": Clade_object}
            root = None

            # open file to parese line iterately
            for i, line in enumerate(pathFile):
                line = line.strip()
                if line.endswith(";"):
                    line = line.rstrip(";")
                line = line.strip().split(";")
                if root is None:
                    root = line[1]
                else:
                    assert root == line[1], "The %d-th line is from a different root"%(i+1)

                # check node iterately, first reverse list, to from leaf to root
                # to make sure every node has a parent node
                leaf2root = line[::-1]

                for j, item in enumerate(leaf2root):
                    # find child_node and parent_node, root node's parent is itself
                    if j == len(line)-1:
                        child_node = item; parent_node=item
                    else:
                        child_node = item; parent_node = leaf2root[j+1]

                    if nodes.has_key(child_node):
                        continue
                    else:
                        # add this node
                        nodes[child_node] = Newick.Clade(name=child_node)
                        # add its parent info
                        nodes[child_node].parent = parent_node

            for node_name, node_clade in nodes.iteritems():
                # find the root node, its parent is itself
                if node_name == node_clade.parent:
                    root_node = node_clade
                    print node_clade
                    print "root node found!! "
                # if node is not root, then find its parent, and add to its parent's clades
                else:
                    parent_node = nodes[node_clade.parent]
                    parent_node.clades.append(node_clade)
                del node_clade.parent

            # transform between output node format
            if node_fmt=="taxid":
                tree = Newick.Tree(root = root_node)
            else:
                assert node_fmt =="sciName", "The node_fmt should be taxid or sciName"
                # convert taxid to sciName
                for node_name, node in nodes.iteritems():
                    node_name = self.get_sciName(node_name)
                    for child in node.clades:
                        if child:
                            child.name = self.get_sciName(child.name)
                root_node.name = self.get_sciName(root_node.name)
                tree = Newick.Tree(root = root_node)

            # write tree to file
            print 'Writing %s tree to %s...' % (out_fmt, outFile)

            bp.write(tree, outFile, out_fmt)


    def taxid2tree(self, taxid_list, out_fmt="newick"):
        """ This function take a list of gi as input, will generate a path for
            for each gi, then construct a newick or phyloxml tree based on these
            gi pathes.

            out_fmt = newick / phyloxml ...
        """
        treeFile = StringIO()

        # get pathes for a list of taxid
        path_list =[";".join([str(item) for item in self.get_path(taxid)])
                    for taxid in taxid_list ]

        # read in pathFile, and store node info into nodes
        nodes = {} # data format {"node_name": Clade_object}
        root = None

        # to parese path iterately
        for i, path in enumerate(path_list):
            line = path.strip().split(";")
            if root is None:
                root = line[0]
            else:
                assert root == line[0], "The %d-th line is from a different root"%(i+1)

            # check node iterately, first reverse list, to from leaf to root
            # to make sure every node has a parent node
            leaf2root = line[::-1]

            for j, item in enumerate(leaf2root):
                # find child_node and parent_node, root node's parent is itself
                if j == len(line)-1:
                    child_node = item; parent_node=item
                else:
                    child_node = item; parent_node = leaf2root[j+1]

                if nodes.has_key(child_node):
                    continue
                else:
                    # add this node
                    nodes[child_node] = Newick.Clade(name=child_node)
                    # add its parent info
                    nodes[child_node].parent = parent_node

        for node_name, node_clade in nodes.iteritems():
            # find the root node, its parent is itself
            if node_name == node_clade.parent:
                root_node = node_clade
                print "root node is %s, constructing tree ..."%(str(node_name))
            # if node is not root, then find its parent, and add to its parent's clades
            else:
                parent_node = nodes[node_clade.parent]
                parent_node.clades.append(node_clade)
            del node_clade.parent


        tree = Newick.Tree(root = root_node)

        bp.write(tree, treeFile, out_fmt)

        treeStr = treeFile.getvalue()

        return treeStr







def main():

    if len(sys.argv) !=2:
        print "\nUsage: python %s[input_centrifuge_output]\n"%sys.argv[0]
        sys.exit(0)

    # with names.dmp and nodes.dmp files
    Tax_ID = TaxIds(names="/data/shengwei/Alteromonas_meta/centrifuge/taxdmp/names.dmp", nodes="/data/shengwei/Alteromonas_meta/centrifuge/taxdmp/nodes.dmp")
    # no names.dmp or nodes.dmp supplied
    #Tax_ID = TaxIds()

    centrifuge_records = CentrifugeRecordParser(sys.argv[1])

    # prepare contigs_dict
    contigs_dict = {}# {contigID:[centrifugeRecord1, centrifugeRecord2, ..]}
    for rec in centrifuge_records:
        if rec.readID in contigs_dict:
            if int(rec.taxID) != 0 and int(rec.score) >= 100:
                contigs_dict[rec.readID].append(rec)
            else:
                continue
        else:
            if int(rec.taxID) != 0 and int(rec.score) >= 100:
                contigs_dict[rec.readID] = [rec]
            else:
                continue

    # write out lca assigned taxonomy for each contig
    with open("contigs_centrifuge_LCA.tsv", "w") as oh:
        oh.write("#ContigID\tLCA_TaxID\tLCA_SciName\tLCA_Rank\tLCA_Path\tMaxScore\tSuperkingdom\tphylum\tclass\torder\tfamily\tgenus\n")
        for contig, hits in contigs_dict.iteritems():

            # get all scores of each hit, find the max score, set 0.98*max as cutoff for LCA
            rec_scores = [int(hit.score) for hit in hits]
            max_score = max(rec_scores)
            cutoff_score = int(max_score * 0.98)
            rec_taxIDs = [int(hit.taxID) for hit in hits if hit.score >= cutoff_score]

            # LCA analysis
            lca_taxid = Tax_ID.get_lca(rec_taxIDs)
            lca_scientificName = Tax_ID.get_sciName(lca_taxid)
            lca_rank = Tax_ID.get_rank(lca_taxid)
            lca_path = Tax_ID.get_path(lca_taxid, toStr=False)

            # get superkingdom
            superkingdom = Tax_ID.get_query_rank_from_path(lca_path, "superkingdom", toStr=True)
            if not superkingdom:
                superkingdom = "None"

            # get phylum
            phylum = Tax_ID.get_query_rank_from_path(lca_path, "phylum", toStr=True)
            if not phylum:
                phylum = "None"

            # get class
            _class = Tax_ID.get_query_rank_from_path(lca_path, "class", toStr=True)
            if not _class:
                _class = Tax_ID.get_query_rank_from_path(lca_path, "subclass", toStr=True)
            if not _class:
                _class = "None"

            # get order
            order = Tax_ID.get_query_rank_from_path(lca_path, "order", toStr=True)
            if not order:
                order = "None"

            # get family
            family = Tax_ID.get_query_rank_from_path(lca_path, "family", toStr=True)
            if not family:
                family = "None"

            # get genus
            genus = Tax_ID.get_query_rank_from_path(lca_path, "genus", toStr=True)
            if not genus:
                genus = "None"



            lca_path = Tax_ID.get_path(lca_taxid, toStr=True)
            oh.write(contig+"\t"+str(lca_taxid)+"\t"+lca_scientificName+"\t"+lca_rank+"\t"+lca_path+"\t"+str(max_score)+"\t"+superkingdom+"\t"+phylum+"\t"+_class+"\t"+order+"\t"+family+"\t"+genus+"\n")



if __name__ == '__main__':
    main()

