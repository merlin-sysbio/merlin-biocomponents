package pt.uminho.ceb.biosystems.merlin.biocomponents.io.writers;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.xml.stream.XMLStreamException;

import org.sbml.jsbml.Annotation;
import org.sbml.jsbml.CVTerm;
import org.sbml.jsbml.CVTerm.Qualifier;
import org.sbml.jsbml.CVTerm.Type;
import org.sbml.jsbml.SBMLDocument;
import org.sbml.jsbml.SBMLException;
import org.sbml.jsbml.xml.XMLNode;
import org.sbml.jsbml.xml.XMLTriple;

import pt.uminho.ceb.biosystems.merlin.biocomponents.container.SBML_Model;
import pt.uminho.ceb.biosystems.merlin.biocomponents.io.SBMLLevelVersion;
import pt.uminho.ceb.biosystems.merlin.biocomponents.io.readers.ContainerBuilder;
import pt.uminho.ceb.biosystems.merlin.database.connector.databaseAPI.CompartmentsAPI;
import pt.uminho.ceb.biosystems.merlin.database.connector.databaseAPI.ModelAPI;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.Connection;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.DatabaseAccess;
import pt.uminho.ceb.biosystems.merlin.utilities.RulesParser;
import pt.uminho.ceb.biosystems.merlin.utilities.Utilities;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.model.CompartmentContainer;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.model.MetaboliteContainer;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.model.ReactionContainer;
import pt.uminho.ceb.biosystems.mew.utilities.datastructures.pair.Pair;


/**
 * @author Oscar
 *
 */
public class SBMLWriter {

	private static final String MERLIN_NAME = "merlin - www.merlin-sysbio.org";
	private static final boolean concatenate = true;
	private boolean isCompartmentalisedModel;
	private Map<String, List<MetaboliteContainer>> reactionMetabolites;
	private Map<String, ReactionContainer> reactions;
	private Map<String, CompartmentContainer> compartments;
	private int compartmentCounter;
	private Map<String, String> compartmentID;
	private String  outsideID;
	private boolean generateFormulae;
	private String biomassEquationID;
	private Map<String,String> metabolites_formula;
	private SBMLLevelVersion levelAndVersion = SBMLLevelVersion.L2V1;
	private DatabaseAccess databaseAccess;
	private String filePath;
	private String sbmlFileID;

	private List<String> genes;
	private SBMLDocument sbmlDocument;

	/**
	 * @param DatabaseAccess
	 * @param filePath
	 * @param sbmlFileID
	 * @param isCompartmentalisedModel
	 * @param generateFormulae
	 * @param biomassEquationID
	 * @param levelAndVersion
	 */
	public SBMLWriter(DatabaseAccess databaseAccess, String filePath, String sbmlFileID, boolean isCompartmentalisedModel,
			boolean generateFormulae, String biomassEquationID, SBMLLevelVersion levelAndVersion) {

		this.databaseAccess = databaseAccess;
		this.compartmentCounter = 1;
		this.compartmentID = new TreeMap<String, String>();
		this.isCompartmentalisedModel = isCompartmentalisedModel;
		this.compartments = new HashMap<String, CompartmentContainer>();
		this.reactionMetabolites = new HashMap<String, List<MetaboliteContainer>>();
		this.reactions = new HashMap<String, ReactionContainer>();
		this.outsideID=null;
		this.generateFormulae = generateFormulae;
		this.metabolites_formula = new TreeMap<String, String>();		
		this.setBiomassEquationName(biomassEquationID);
		this.filePath = filePath;
		this.sbmlFileID = sbmlFileID;
		this.levelAndVersion = levelAndVersion;

		this.genes = new ArrayList<>();
	}


	/**
	 * 
	 */
	public void getDataFromDatabase() {

		String aux="";

		if(this.isCompartmentalisedModel)
			aux = aux.concat(" NOT originalReaction");
		else
			aux = aux.concat(" originalReaction");

		this.getReactions(databaseAccess, aux);
		this.getStoichiometry(databaseAccess, aux);
	}

	/**
	 * <p>This method converts the model data to the SBML</p>
	 * <p>native format and returns it as an SBML <code>String</code></p> 
	 * @param addAllNotes 
	 * 
	 * @return <code>String</code> representation of the SBML model.
	 * @throws Exception 
	 */
	public void toSBML(boolean addAllNotes) throws Exception{

		SBML_Model sbmlModel = new SBML_Model(sbmlFileID.replaceAll("[^a-zA-Z0-9]", ""), levelAndVersion );
		this.getCompartments(databaseAccess, sbmlModel);
		this.buildModel(sbmlModel, addAllNotes,databaseAccess);

		this.sbmlDocument = this.createSBMLDocument(sbmlModel, filePath);

		if(this.generateFormulae)
			this.generateFormulae(filePath);
	}


	/**
	 * @param sbmlModel
	 * @param filePath
	 * @throws FileNotFoundException
	 * @throws SBMLException
	 * @throws XMLStreamException
	 */
	public SBMLDocument createSBMLDocument(SBML_Model sbmlModel, String filePath) throws FileNotFoundException, SBMLException, XMLStreamException {

		SBMLDocument document = new SBMLDocument(levelAndVersion.getLevel(), levelAndVersion.getVersion());
		document.setModel(sbmlModel.getModel());
		document.addDeclaredNamespace("xmlns:html", "http://www.w3.org/1999/xhtml");

		org.sbml.jsbml.SBMLWriter sbmlwriter = new org.sbml.jsbml.SBMLWriter();
		sbmlwriter.setProgramName(MERLIN_NAME);
		String merlinVersion = Utilities.getMerlinVersion();

		if(merlinVersion!=null)
			sbmlwriter.setProgramVersion(merlinVersion);

		OutputStream out = new FileOutputStream(filePath);
		sbmlwriter.write(document, out);
		try {
			out.close();
		} 
		catch (IOException e) {

			e.printStackTrace();
		}
		return document;
	}

	/**
	 * @param link
	 * @param conditions
	 */
	private void getReactions(DatabaseAccess link, String conditions) {

		Statement stmt;
		try {

			Connection connection = new Connection(link);

			stmt = connection.createStatement();

			Map<Integer, Pair<String, String>> genesData = ModelAPI.getGenesFromDatabase(stmt);

			Map<String,ArrayList<String>> result = ModelAPI.getReactions(stmt, conditions);

			for(String idReaction : result.keySet()){

				ArrayList<String> reactionData = result.get(idReaction);

				ReactionContainer reactionContainer = null;

				if(this.reactions.containsKey(idReaction)) {

					System.out.println("same reaction "+ idReaction);
				}
				else {

					reactionContainer = new ReactionContainer(idReaction, reactionData.get(0), reactionData.get(1), Boolean.valueOf(reactionData.get(2)), reactionData.get(3), reactionData.get(0));
					if(reactionData.get(5)!= null && !reactionData.get(5).equalsIgnoreCase("null") && !reactionData.get(5).isEmpty())
						reactionContainer.setLowerBound(Double.parseDouble(reactionData.get(5)));
					if(reactionData.get(6)!= null && !reactionData.get(6).equalsIgnoreCase("null")&& !reactionData.get(6).isEmpty())
						reactionContainer.setUpperBound(Double.parseDouble(reactionData.get(6)));

					if(reactionData.get(4)!=null)
						reactionContainer.setNotes(reactionData.get(4));

					if(reactionData.get(7)!=null && !reactionData.get(7).isEmpty()){
						String geneRules = RulesParser.getGeneRule(reactionData.get(7), genesData);
						reactionContainer.setGeneRule(geneRules);
					}
				}

				reactionContainer = this.reactions.put(idReaction, reactionContainer);
			}

			ArrayList<String[]> result2 = ModelAPI.getReactionHasEnzymeData3(stmt);

			for(int i=0; i<result2.size(); i++){
				String[] list = result2.get(i);

				if(this.reactions.containsKey(list[0]) && !list[2].contains(".-")) {

					ReactionContainer reactionContainer = this.reactions.get(list[0]);

					Set<String> enzymeSet = new TreeSet<String>();
					if(reactionContainer.getEnzymes()!=null)
						enzymeSet = reactionContainer.getEnzymes();

					if(list[2]!=null) {

						enzymeSet.add(list[2]);
						reactionContainer.setEnzymes(enzymeSet);
					}
					reactionContainer = this.reactions.put(list[0], reactionContainer);
				}
			}

			result2 = ModelAPI.getPathwayHasReactionData(stmt);

			for(int i=0; i<result2.size(); i++) {

				String[] list = result2.get(i);

				if(this.reactions.containsKey(list[0])) {

					ReactionContainer reactionContainer = this.reactions.get(list[0]);

					Set<String> pathways = new TreeSet<String>();
					if(reactionContainer.getPathways()!=null) {

						pathways = reactionContainer.getPathways();
					}
					if(list[1]!=null) {

						pathways.add(list[2]);
						reactionContainer.setPathways(pathways);
					}
					reactionContainer = this.reactions.put(list[0], reactionContainer);
				}
			}

			result2 = ModelAPI.getReactionGenes(stmt);

			for(int i=0; i<result2.size(); i++){
				String[] list = result2.get(i);

				if(this.reactions.containsKey(list[0]) && !list[3].contains(".-")) {

					ReactionContainer reactionContainer = this.reactions.get(list[0]);

					String locus= list[2], geneName = null;

					if(list[1]!=null)
						geneName = list[1].replace(" ","").replace(",","_").replace("/","_").replace("\\","_").trim();//replace("-","_").trim();

					//					this.addGeneCI(locus, geneName);
					reactionContainer.addGene(locus, geneName);
					reactionContainer = this.reactions.put(list[0], reactionContainer);

				}
			}
			stmt.close();
		}
		catch (SQLException e) {

			e.printStackTrace();
		}
	}

	/**
	 * @param link
	 * @param conditions
	 */
	private void getStoichiometry(DatabaseAccess link, String conditions) {

		Statement stmt;
		try {

			Connection connection = new Connection(link);

			stmt = connection.createStatement();

			ArrayList<String[]> result = ModelAPI.getStoichiometry(conditions, stmt);

			for(int i=0; i<result.size(); i++){
				String[] list = result.get(i);

				if(this.reactions.containsKey(list[1])) {

					if(!list[4].contains("m") && !list[4].contains("n")) {

						List<MetaboliteContainer> metabolitesContainer = new ArrayList<MetaboliteContainer>();

						if(this.reactionMetabolites.containsKey(list[1]))
							metabolitesContainer = this.reactionMetabolites.get(list[1]);

						MetaboliteContainer metabolite = new MetaboliteContainer(Integer.parseInt(list[2]), list[6], list[7], list[4], list[5], list[3]);
						metabolite.setEntryID(list[8]);

						metabolitesContainer.add(metabolite);

						this.reactionMetabolites.put(list[1], metabolitesContainer);
					}
					else {

						this.reactionMetabolites.remove(list[1]);
						this.reactions.remove(list[1]);
					}
				}
			}
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}


	/**
	 * @param link
	 * @param conditions
	 */
	private void getCompartments(DatabaseAccess link, SBML_Model sbml) {

		Statement stmt;
		try {

			Connection connection = new Connection(link);

			stmt = connection.createStatement();

			ArrayList<String[]> result = CompartmentsAPI.getCompartments2(stmt);

			for(int i=0; i<result.size(); i++){
				String[] list = result.get(i);

				CompartmentContainer compartmentContainer = new CompartmentContainer(list[0], list[1],  list[2]);
				this.compartments.put(list[0], compartmentContainer);

				if(( list[1].equalsIgnoreCase("extracellular") && this.isCompartmentalisedModel) || ( list[1].equalsIgnoreCase("outside") && !this.isCompartmentalisedModel))
					this.outsideID=this.getCompartmentID( list[0], sbml);
			}

		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

	/**
	 * @param compartment
	 * @return
	 */
	private String getCompartmentID(String compartment, SBML_Model sbml) {

		if(!this.compartmentID.containsKey(compartment)) {

			String id = SBMLWriter.buildID("C_", this.compartmentCounter);
			sbml.addCompartment(id, this.compartments.get(compartment).getName(), this.outsideID);
			this.compartmentID.put(compartment, id);
			this.compartmentCounter++ ;

		}
		return this.compartmentID.get(compartment);
	}

	/**
	 * @param addAllNotes 
	 * @throws SQLException 
	 * 
	 */
	public void buildModel(SBML_Model sbmlModel, boolean addAllNotes, DatabaseAccess link) throws SQLException {

		int reactionsCounter = 1 ;
		int metabolitesCounter = 1;
		int enzymesCounter = 1;

		Map<String,String> enzymesID = new TreeMap<String, String>();
		Map<String,String> compoundCompartmentID = new TreeMap<String, String>();

		for(String reaction_id : this.reactions.keySet()) {

			ReactionContainer reaction = this.reactions.get(reaction_id);
			String name = reaction.getName()+"__("+reaction.getEquation().replace(" ", "")+")";

			String rid = SBMLWriter.buildID("R_", reactionsCounter)/*+"_"+reaction.getName().replace(" ", "_").replace("\t", "_").replace("-", "_")*/;

			boolean biomassEquation = reaction.getName().equalsIgnoreCase(this.getBiomassEquationID());

			double upper_bound = 999999;
			double lower_bound = 0;

			if(reaction.getLowerBound()!= null)
				lower_bound = reaction.getLowerBound();
			else 
				if(reaction.isReversible())
					lower_bound = -999999;

			if(reaction.getUpperBound()!= null)
				upper_bound = reaction.getUpperBound();


			if(reaction.getGenes()!= null)
				for(Pair<String,String> gene : reaction.getGenes())			
					if(!genes.contains(gene.getA()))
						genes.add(gene.getA());

			String geneRule = reaction.getGeneRule();

			//			System.out.println(reaction.getGeneRule());

			//			if(reaction.getGeneRule()!=null && !reaction.getGeneRule().isEmpty()){
			//				List<List<Pair<String,String>>> rule = ModelAPI.parseBooleanRule(reaction.getGeneRule(), stmt);
			//				geneRule = Utilities.parseRuleListToString(rule);
			//			}

			if(geneRule!= null && !geneRule.isEmpty())
				reaction.setGeneRule(geneRule);


			XMLNode reactionNote = SBMLWriter.getReactionNote(reaction.getGenes(), reaction.getEnzymes(),
					reaction.getPathways(), null, reaction.getNotes(), addAllNotes, reaction.getGeneRule(), this.levelAndVersion);


			sbmlModel.addReaction(rid, name, reaction.getName(), reaction.isReversible(), reactionNote, lower_bound, upper_bound, biomassEquation);
			reactionsCounter++ ;

			if(this.reactions.get(reaction_id).getEnzymes() != null) {

				for(String enzyme : this.reactions.get(reaction_id).getEnzymes()) {

					String enzyme_surrogate = enzyme.concat("_").concat(this.getCompartmentID(reaction.getLocalisation(), sbmlModel));
					String eid;
					if(enzymesID.containsKey(enzyme_surrogate)) {

						eid = enzymesID.get(enzyme_surrogate);
					}
					else {

						eid = SBMLWriter.buildID("E_", enzymesCounter);

						sbmlModel.addSpecies(eid, enzyme, this.getCompartmentID(reaction.getLocalisation(), sbmlModel), SBMLWriter.getAnnotation(enzyme, ""));
						enzymesID.put(enzyme_surrogate, eid);
						enzymesCounter++;
					}
					sbmlModel.setReactionEnzyme(eid, rid);
				}
			}

			if(this.reactionMetabolites.get(reaction_id) != null) {

				for(MetaboliteContainer metabolite : this.reactionMetabolites.get(reaction_id)) {

					String metabolite_surrogate = metabolite.getMetaboliteID()+"".concat("_").concat(metabolite.getCompartment_name());
					String mid;
					if(compoundCompartmentID.containsKey(metabolite_surrogate)) {

						mid = compoundCompartmentID.get(metabolite_surrogate);
					}
					else {

						mid = SBMLWriter.buildID("M_", metabolitesCounter);
						compoundCompartmentID.put(metabolite_surrogate,mid);

						String sbmlName="";

						if(metabolite.getEntryID() != null) {

							sbmlName=sbmlName.concat(metabolite.getEntryID());
						}

						if(metabolite.getName() != null) {

							if(sbmlName!=null && !sbmlName.isEmpty())
								sbmlName=sbmlName.concat("_");

							sbmlName=sbmlName.concat(metabolite.getName());
						}

						if(metabolite.getFormula() != null) {

							if(sbmlName!=null && !sbmlName.isEmpty())
								sbmlName=sbmlName.concat("_");

							sbmlName=sbmlName.concat(metabolite.getFormula());
						}
						sbmlModel.addSpecies( mid, sbmlName, this.getCompartmentID(metabolite.getCompartment_name(), sbmlModel), SBMLWriter.getAnnotation(metabolite.getEntryID(), metabolite.getFormula()));

						metabolitesCounter++ ;

						metabolites_formula.put(mid, metabolite.getFormula());
					}

					sbmlModel.setReactionCompound(mid, rid, new Double(metabolite.getStoichiometric_coefficient()));

				}
			}
			else {

				System.err.println(reaction_id +"\t" +reaction.getEntryID());
			}
		}

	}

	/**
	 * @param filePath
	 */
	public void generateFormulae(String filePath){

		if(this.metabolites_formula.size()>0) {

			try {

				FileWriter fstream = new FileWriter(filePath.replace(".xml","")+"_formulae.txt");
				BufferedWriter out = new BufferedWriter(fstream);

				for(String mid:metabolites_formula.keySet()) {

					out.write(mid+"\t"+metabolites_formula.get(mid)+"\n");
				}
				out.close();
			}
			catch (Exception e){

				System.err.println("Error: " + e.getMessage());
			}
		}
	}

	/**
	 * @param type
	 * @param counter
	 * @return
	 */
	private static String buildID(String type, int counter) {

		if(counter<10 ){return type.concat("0000"+counter);}
		if(counter<100) {return type.concat("000"+counter);}
		if(counter<1000) {return type.concat("00"+counter);}
		if(counter<10000) {return type.concat("0"+counter);}

		return type.concat("_"+counter);
	}

	/**
	 * @param urn_id
	 * @param formula
	 * @return
	 */
	public static Annotation getAnnotation(String urn_id, String formula) {

		Annotation annotation = new Annotation();
		CVTerm cvTerm = new CVTerm(Type.BIOLOGICAL_QUALIFIER, Qualifier.BQB_IS);
		annotation.addCVTerm(cvTerm);

		if(urn_id!=null) {

			if(urn_id.startsWith("C")) {

				if(formula!=null) {

					cvTerm.addResource("FORMULA:" + formula.toUpperCase());
				}
				cvTerm.addResource("urn:miriam:kegg.compound:" + urn_id);
				//cvTerm.addResource("<url:element>http://www.genome.jp/dbget-bin/www_bget?gl:"+urn_id+ "</url:element>");
				//annotation.appendNoRDFAnnotation("http://www.genome.jp/dbget-bin/www_bget?gl:"+urn_id);
				cvTerm.addResource("http://www.genome.jp/dbget-bin/www_bget?cpd:"+urn_id);
			}

			else if(urn_id.startsWith("G")) {

				if(formula!=null) {

					cvTerm.addResource("FORMULA: " + formula.toUpperCase());
				}
				cvTerm.addResource("urn:miriam:kegg.glycan:" + urn_id);
				cvTerm.addResource("http://www.genome.jp/dbget-bin/www_bget?gl:"+urn_id);
			}

			else if(urn_id.startsWith("D")) {

				if(formula!=null) {

					cvTerm.addResource("FORMULA: " + formula.toUpperCase());
				}
				cvTerm.addResource("urn:miriam:kegg.drugs:" + urn_id);
				cvTerm.addResource("http://www.genome.jp/dbget-bin/www_bget?dr:"+urn_id);
			}

			else if (urn_id.contains("#")) {

				cvTerm.addResource("urn:miriam:tcdb:" + urn_id);
				cvTerm.addResource("http://www.tcdb.org/search/result.php?tc="+urn_id);
			}

			else if (urn_id.contains(".")) {

				cvTerm.addResource("urn:miriam:ec-code:" + urn_id);
				cvTerm.addResource("http://www.genome.jp/dbget-bin/www_bget?ec:"+urn_id);
			}
		}
		return annotation;

	}

	/**
	 * @param genes
	 * @param proteins
	 * @param pathways
	 * @param proteinClass
	 * @param notes
	 * @param addAllNotes
	 * @param geneRule
	 * @param levelAndVersion
	 * @return
	 */
	public static XMLNode getReactionNote(Set<Pair<String,String>> genes, Set<String> proteins, Set<String> pathways, Set<String> proteinClass, String notes, boolean addAllNotes, String geneRule, SBMLLevelVersion levelAndVersion) {

		XMLNode note = new XMLNode(new XMLTriple("notes"));
		note.addNamespace("html");

		String genesNotes = RulesParser.processReactionGenes(genes, concatenate);
		String proteinsNotes = ContainerBuilder.processReactionProteins(proteins);
		String pathwaysNotes = ContainerBuilder.processReactionPathways(pathways);
		String enzymesNotes = ContainerBuilder.processReactionProteinClass(proteinClass);

		Set<String> notesList = ContainerBuilder.processReactionNotes(notes);

		String newGeneRule = geneRule;
		if(newGeneRule!= null && !levelAndVersion.equals(SBMLLevelVersion.L3V1))			
			newGeneRule = newGeneRule.replaceAll(" \\(", "_").replaceAll("\\)", "");

		//		System.out.println(newGeneRule);

		Set<XMLNode> gener = ContainerBuilder.getGeneRules(genesNotes, notesList, addAllNotes, newGeneRule);
		Set<XMLNode> proteinr = ContainerBuilder.getProteinRules(proteinsNotes, notesList, addAllNotes);
		Set<XMLNode> proteinc = ContainerBuilder.getPathwaysRules(pathwaysNotes, notesList, addAllNotes);
		Set<XMLNode> sub = ContainerBuilder.getProteinClassRules(enzymesNotes, notesList, addAllNotes);

		for(XMLNode x : gener)
			note.addChild(x);

		for(XMLNode x : proteinr)
			note.addChild(x);

		for(XMLNode x : proteinc)
			note.addChild(x);

		for(XMLNode x : sub)
			note.addChild(x);

		return note;		
	}


	/**
	 * @param notes
	 * @return
	 */
	public static Set<String> processReactionNotes(String notes) {

		Set<String> reactionNotes = new HashSet<>();
		//<html:p>GENE_ASSOCIATION: </html:p>
		//<html:p>PROTEIN_ASSOCIATION: "</html:p>
		//<html:p>SUBSYSTEM: S_Tyrosine__Tryptophan__and_Phenylalanine_Metabolism</html:p>
		//<html:p>PROTEIN_CLASS: 1.13.11.27</html:p>
		//<html:p>GENE_ASSOCIATION: ( YOL096C  and  YDR204W  and  YML110C  and  YGR255C  and  YOR125C  and  YGL119W  and  YLR201C )</html:p>
		//<html:p>PROTEIN_ASSOCIATION: ( Coq3-m and Coq4-m and Coq5-m and Coq6-m and Coq7-m and Coq8-m and Coq9-m )"</html:p>
		//<html:p>SUBSYSTEM: S_Quinone_Biosynthesis</html:p>
		//<html:p>PROTEIN_CLASS: </html:p>


		if(notes!=null && !notes.trim().isEmpty())
			if(notes.contains("|"))
				for(String note : notes.split(" \\| "))
					reactionNotes.add(note.replace(",","_").replace(" ","_").replace(":_",": ").replace("_AND_"," and ").replace("_OR_"," or ").trim());	 //.replace(")","").replace("(","_")

		return reactionNotes; 
	}


	/**
	 * @return the biomassEquationID
	 */
	public String getBiomassEquationID() {
		return biomassEquationID;
	}

	/**
	 * @param biomassEquationID the biomassEquationID to set
	 */
	public void setBiomassEquationName(String biomassEquationID) {
		this.biomassEquationID = biomassEquationID;
	}


	/**
	 * @return
	 */
	public SBMLDocument getDocument() {

		return this.sbmlDocument;
	}

}
