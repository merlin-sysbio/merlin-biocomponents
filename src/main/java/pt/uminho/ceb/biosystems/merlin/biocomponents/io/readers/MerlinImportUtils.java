package pt.uminho.ceb.biosystems.merlin.biocomponents.io.readers;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;

import pt.uminho.ceb.biosystems.merlin.biocomponents.io.Enumerators.ModelSources;
import pt.uminho.ceb.biosystems.merlin.utilities.RulesParser;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.model.CompartmentContainer;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.model.EnzymeContainer;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.model.GeneContainer;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.model.MetaboliteContainer;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.model.PathwaysHierarchyContainer;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.model.ReactionContainer;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.modelSeed.ModelSeedCompoundsDB;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.modelSeed.ModelSeedPathwaysDB;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.modelSeed.ModelSeedReactionsDB;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.CompartmentCI;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.GeneCI;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.MetaboliteCI;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.ReactionCI;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.ReactionTypeEnum;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.StoichiometryValueCI;
import pt.uminho.ceb.biosystems.mew.utilities.grammar.syntaxtree.AbstractSyntaxTreeNode;
import pt.uminho.ceb.biosystems.mew.utilities.math.language.mathboolean.DataTypeEnum;
import pt.uminho.ceb.biosystems.mew.utilities.math.language.mathboolean.IValue;


public class MerlinImportUtils {

	private MerlinSBMLContainer cont;
	private ConcurrentLinkedQueue<MetaboliteContainer> resultMetabolites;
	private ConcurrentLinkedQueue<EnzymeContainer> resultEnzymes;
	private ConcurrentLinkedQueue<ReactionContainer> resultReactions;
	private ConcurrentLinkedQueue<GeneContainer> resultGenes;
	private ConcurrentLinkedQueue<CompartmentContainer> resultCompartments;
	private Map<String, String> enzymeNames;
	private ConcurrentLinkedQueue<PathwaysHierarchyContainer> resultPathwaysHierarchy;
	private ModelSeedReactionsDB reactionsData;
	private ModelSeedPathwaysDB keggPathwaysData;
	private ModelSeedCompoundsDB metabolitesData;
	private Map<String, String> metaboliteCompartments;
	private List<String> transportReactions;
	private List<String> drains;
	private ModelSources modelSource;
	private Set<String> metabolitesSet;
//	private Map<String,Integer> genesIds;
	private Map<String,String> compartments;
	private Map<String,Set<String>> reactionEnzymes;
//	private boolean addSpontaneousReactions;


//	public MerlinImportUtils(MerlinSBMLContainer sbml3Container, ModelSources source, String level) throws SQLException {
//		
//		this(sbml3Container, source);
//	}

	public MerlinImportUtils(MerlinSBMLContainer container, ModelSources source) throws SQLException{

		this.cont = container;
		this.modelSource = source;
		this.resultReactions = new ConcurrentLinkedQueue<>();
		this.resultGenes = new ConcurrentLinkedQueue<>();
		this.resultMetabolites = new ConcurrentLinkedQueue<>();
		this.resultEnzymes = new ConcurrentLinkedQueue<>();
		this.resultCompartments = new ConcurrentLinkedQueue<>();
		this.enzymeNames = new HashMap<>();
		this.reactionEnzymes = new HashMap<>();
		
		if(this.modelSource.equals(ModelSources.MODEL_SEED)){
			this.metabolitesData = new ModelSeedCompoundsDB();
			this.reactionsData = new ModelSeedReactionsDB();
		}	
		
		this.keggPathwaysData = new ModelSeedPathwaysDB();
		this.metaboliteCompartments = new HashMap<>();
		this.metabolitesSet = new HashSet<>();
		this.compartments = new HashMap<>();
		
//		if(statement!=null)
//			this.genesIds = ModelAPI.getGeneIds(statement);
		
		this.resultPathwaysHierarchy = new ConcurrentLinkedQueue<>();
		this.transportReactions = new ArrayList<>(cont.getReactionsByType(ReactionTypeEnum.Transport));
		this.drains = new ArrayList<>();
//		this.addSpontaneousReactions = false;

		
		System.out.println("Reading compartments...");
		readCompartments();
		System.out.println("Reading metabolites...");
		readMetabolites();
		System.out.println("Reading pathways...");
		readPathways();
		System.out.println("Reading genes...");
		readGenes();
		System.out.println("Reading reactions...");
		readReactions();
		System.out.println("Reading enzymes...");
		readEnzymes();
		
	}
	
	

	/**
	 * 
	 */
	private void readGenes(){

		//		GENE CONTAINER
		//		private String entry;	X
		//		private String name;	X
		//		private	List<String> dblinks;
		//		private List<String> orthologues;
		//		private	List<String> genes;
		//		private	List<String> modules;
		//		private String chromosome_name;
		//		private String  left_end_position, right_end_position, aasequence, aalength, ntsequence, ntlength;

		Map<String, GeneCI> genes = cont.getGenes();
		
		if(!genes.isEmpty()){

			for(GeneCI gene : genes.values()){

				String geneID = gene.getGeneId();
				
//				System.out.println("GeneID: "+geneID);
				
				if(geneID!=null && !geneID.isEmpty()){
				 
					GeneContainer geneContainer = new GeneContainer(geneID);
					if(gene.getGeneName()!=null)
						geneContainer.setName(gene.getGeneName());
					
					this.resultGenes.add(geneContainer);
				}
			}
		}
	}


	/**
	 * 
	 */
	private void readMetabolites(){

		//		METABOLITE CONTAINER
		//		private int metaboliteID;
		//		private String name;	X
		//		private String formula;	X
		//		private String stoichiometric_coefficient;
		//		private String numberofchains;
		//		private String compartment_name;	X
		//		private String entryID;
		//		private String molecular_weight;	X
		//		private List<String> names;
		//		private List<String> enzymes;
		//		private List<String> reactions;	X
		//		private Map<String, String> pathways;
		//		private List<String> dblinks;
		//		private List<String> same_as;

		Map<String, MetaboliteCI> metabolites = cont.getMetabolites();
		
		
//		Set<String> coreMetabolites = new HashSet<>();
//		if(this.metabolitesData!=null)
//			coreMetabolites = this.metabolitesData.getCoreCompounds();
		
		for(String metID : metabolites.keySet()){
			
			if(!metID.endsWith("_b")){

				MetaboliteCI metabolite = metabolites.get(metID);
//				MetaboliteContainer metContainer = new MetaboliteContainer(metID);
				
				String metaboliteID = processSBMLMetaboliteID(metID);//.split("_")[1];
				
				String metaboliteCompartmentID = cont.getMetaboliteCompartments(metID).toArray()[0].toString();
				String compartmentName = this.compartments.get(metaboliteCompartmentID);
//				if(metaboliteCompartmentID!=null && !metaboliteCompartmentID.isEmpty())
//					compartmentName = this.compartments.get(metaboliteCompartmentID);
				
				this.metaboliteCompartments.put(metID, compartmentName);


				if(!this.metabolitesSet.contains(metaboliteID)){
					
					MetaboliteContainer metContainer = new MetaboliteContainer(metaboliteID);

					//			String[] nameElems = metabolite.getName().split("_");
					//			String name = "";

					//			for(String elem : nameElems){
					//
					//				try {
					//
					//					int prefix = Integer.parseInt(elem);
					//
					//					name = name.concat(prefix+"").concat(",");
					//
					//				} 
					//				catch (NumberFormatException e) {
					//
					//					name = name.concat(set[i].trim());
					//
					//					metabolites.add(compound);
					//
					//					compound = "";
					//				}
					//			}
					
					//METABOLITE COMPARTMENT
					metContainer.setCompartment_name(compartmentName);

					//METABOLITE NAME, FORMULA AND MOLECULAR WEIGHT
					if(this.metabolitesData!=null && this.metabolitesData.existsCompoundID(metaboliteID)){
						
						metContainer.setName(this.metabolitesData.getCompoundName(metaboliteID));
						metContainer.setFormula(this.metabolitesData.getCompoundFormula(metaboliteID));
						metContainer.setMolecular_weight(this.metabolitesData.getCompoundMolecularWeight(metaboliteID));
						
//						coreMetabolites.remove(metaboliteID);
					}
					else{
						metContainer.setName(metabolite.getName());
						metContainer.setFormula(metabolite.getFormula());
						Double mass = metabolite.getMass();
						if(mass != null && mass!=0)
							metContainer.setMolecular_weight(mass.toString());

					}

					//METABOLITE REACTIONS
					Set<String> metaboliteReactions = metabolite.getReactionsId();
					List<String> reactionsList = new ArrayList<>(metaboliteReactions);
					metContainer.setReactions(reactionsList);

					//METABOLITE SYMNONYMS
					if(metabolite.getSymnonyms() != null)
						metContainer.setSame_as(metabolite.getSymnonyms());
					
					this.metabolitesSet.add(metaboliteID);
					this.resultMetabolites.add(metContainer);
				}
			}
			else{
				Set<String> reactionsList = cont.getMetabolite(metID).getReactionsId();
				
				for(String drainID : reactionsList)
					if(!drains.contains(drainID))
						drains.add(drainID);
			}
		}
		
//		for(String coreMetabolite : coreMetabolites){
//			
//			MetaboliteContainer metaboliteCont = new MetaboliteContainer(coreMetabolite);
//			
//			metaboliteCont.setName(this.metabolitesData.getCompoundName(coreMetabolite));
//			metaboliteCont.setFormula(this.metabolitesData.getCompoundFormula(coreMetabolite));
//			metaboliteCont.setMolecular_weight(this.metabolitesData.getCompoundMolecularWeight(coreMetabolite));
//			
//			this.resultMetabolites.add(metaboliteCont);
//		}
	}

	/**
	 * 
	 */
	private void readReactions(){

		//		REACTION CONTAINER
		//		private String entryID;		X
		//		private String reactionID;	X
		//		private boolean reversible, inModel;	X	X
		//		private Double lowerBound, upperBound;	X	X
		//		private String name, localisation, notes;	X	X	X
		//		private List<String> names;
		//		private List<String> dblinks;
		//		private String equation;	X
		//		private	Map<String, String[]> reactantsStoichiometry;	X
		//		private	Map<String, String[]> productsStoichiometry;	X
		//		private Set<String> enzymes, comments, genes, pathways;   X  -  X    X
		//		private	Map<String, String> pathwaysMap; 	-
		//		private String geneRule;  X

		Map<String, ReactionCI> reactions = cont.getReactions();
		

		for(String idReaction : reactions.keySet()){
			
			String reactionID = idReaction;
			
			if(idReaction.startsWith("R_"))
				reactionID = idReaction.substring(2);
				
			if(reactionID.contains("_")){
				
				reactionID = reactionID.substring(0, reactionID.lastIndexOf("_"));
				
				String[] splitedID = reactionID.split("_");

//				String id = "";
				if(splitedID.length>0) {
					for(int i=0; i<splitedID.length; i++){
						
						if(modelSource.equals(ModelSources.MODEL_SEED)) 
							if(splitedID[i].matches("rxn\\d{5}"))
								reactionID = splitedID[i];
						
	//					else if(splitedID[i].matches("R\\d{5}"))
	//						id = splitedID[i];
	//					
	//					else
	//						id.concat("_").concat(splitedID[i]);
					}
				}
			}
			
			ReactionContainer reactionContainer = new ReactionContainer(reactionID);
			ReactionCI reaction = reactions.get(idReaction);
			
			Map <String, StoichiometryValueCI> reactants = reaction.getReactants();
			Map <String, StoichiometryValueCI> products = reaction.getProducts();
			
			reactionContainer.setReversible(reaction.isReversible());
			
			String externalID = "";
			
			if(this.reactionsData!=null && this.reactionsData.existsReactionID(reactionID)){
				
				if(this.reactionsData.getReactionAbbreviation(reactionID).contains("_"))
					externalID = this.reactionsData.getReactionAbbreviation(reactionID).split("_")[1];
				else
					externalID = this.reactionsData.getReactionAbbreviation(reactionID);
				
				if(externalID.matches("R\\d{5}")){
					
//					try {
//						Integer.parseInt(externalID.substring(1));	//verify we have an ID number with KeggID format ("R"+ numbers)
						reactionContainer.setEntryID(externalID);
//					} 
//					catch (NumberFormatException e) {
//					}
				}
				
				if(reaction.isReversible() == reactionsData.isReactionReversible(reactionID)){
					
//					String splitChar = "";
//					if(!reactionsData.isReactionReversible(reactionID))
//						splitChar = reactionsData.getReactionDirection(reactionID);
//					else
//						splitChar = "<=>";
					
					String equation = this.reactionsData.getReactionEquation(reactionID);//.replaceAll("[\\[\\(]\\s*\\d*\\s*[\\]\\)]", "").trim();
					
//					if(equation.trim().split(splitChar).length == 2)
//						reactionContainer.setEquation(equation);
					if(!equation.isEmpty())
						reactionContainer.setEquation(equation);
				}
			}
			
			//EQUATION
			if(reactionContainer.getEquation() == null) {
				
				String direction = null;

				if(reaction.isReversible())
					direction = " <=> ";
				else
					direction = " => ";

				String equation = "";
				Map<String, MetaboliteCI> metabolites = cont.getMetabolites();

				if(!reactants.isEmpty() && !products.isEmpty()){
					for(String reactantID : reactants.keySet()){

						String reactant = metabolites.get(reactantID).getName();
						equation = equation.concat(reactant).concat(" + ");
					}

					equation = equation.substring(0, equation.lastIndexOf(" + "));
					equation = equation.concat(direction);

					for(String productID : products.keySet()){

						String product = metabolites.get(productID).getName();
						equation = equation.concat(product).concat(" + ");
					}

					equation = equation.substring(0, equation.lastIndexOf(" + "));
					equation.replaceAll("[\\[\\(]\\s*\\d*\\s*[\\]\\)]", "").trim();
					reactionContainer.setEquation(equation);
				}
				else{
					reactionContainer.setEquation(reaction.getName());
				}
			}
			
			reactionContainer.setReactionID(reactionID);
			for(String geneID : reaction.getGenesIDs())
				reactionContainer.addGene(geneID, cont.getGene(geneID).getGeneName());
			reactionContainer.setEnzymes(reaction.getProteinIds());
			reactionContainer.setInModel(true);
			
			//REACTION COMPARTMENT
			if(!reaction.identifyCompartments().isEmpty())
				reactionContainer.setLocalisation(cont.getCompartment(reaction.identifyCompartments().toArray()[0].toString()).getName().split("_")[0]);
			
//			GENES RULES
			if(reaction.getGeneRule()!=null && reaction.getGeneRule().size()!=0){
				
				AbstractSyntaxTreeNode<DataTypeEnum, IValue> treeRoot = reaction.getGeneRule().getRootNode();
				List<String> geneCombinations = RulesParser.getGeneRuleCombinations(treeRoot);
				
//				String geneRule = RulesParser.getGeneRuleString(geneCombinations, this.genesIds);
				String geneRule = RulesParser.getOR_geneRulesList2String(geneCombinations);
				
				reactionContainer.setGeneRule(geneRule);
			}
			
			//ENZYMES
			Set<String> enzymes = reaction.getEcNumbers();
			reactionContainer.setEnzymes(enzymes);
			String name = reaction.getName();
			reactionContainer.setName(name);
			this.enzymeNames.put(reactionID,name);
			this.reactionEnzymes.put(reactionID,enzymes);

			//BOUNDS
			if(cont.getDefaultEC().containsKey(idReaction)){
				reactionContainer.setLowerBound(cont.getDefaultEC().get(idReaction).getLowerLimit());
				reactionContainer.setUpperBound(cont.getDefaultEC().get(idReaction).getUpperLimit());
			}
			
			//REACTANTS STOICHIOMETRY
			Map<String, String[]> reactantsStoichiometry = new HashMap<>();

			for(String reactantID : reactants.keySet()){
				
				if(!reactantID.endsWith("_b")){
				
	//				String kBaseReactantID = reactantID.split("_")[1];
					String[] stoichiometryValue = new String[3];
					stoichiometryValue[0] = Double.toString(-reactants.get(reactantID).getStoichiometryValue());
					stoichiometryValue[1] = "";
					stoichiometryValue[2] = this.metaboliteCompartments.get(reactantID);
					
					//verify if compartments of metabolites in reaction matches
					String compID = null;
					for(String component : reactantID.split("_"))
						if(this.compartments.keySet().contains(component))
							compID=component;
					
					if(compID!=null && !compID.isEmpty()){
						
						String compartmentName = this.compartments.get(compID);
//						if(compartmentName.contains("_"))
//							compartmentName = compartmentName.split("_")[0];
						
						if(!compartmentName.equals(this.metaboliteCompartments.get(reactantID))){
							if(reactionContainer.getNotes()==null || reactionContainer.getNotes().isEmpty())
								reactionContainer.setNotes("verify the reactants compartments for this reaction");
						}
					}
						
					reactantsStoichiometry.put(reactantID, stoichiometryValue);
				}
			}

			reactionContainer.setReactantsStoichiometry(reactantsStoichiometry);

			//PRODUCTS STOICHIOMETRY
			Map<String, String[]> productsStoichiometry = new HashMap<>();

			for(String productID : products.keySet()){
				
				if(!productID.endsWith("_b")){
				
	//				String kBaseProductID = productID.split("_")[1];
					String[] stoichiometryValue = new String[3];
					stoichiometryValue[0] = Double.toString(products.get(productID).getStoichiometryValue());
					stoichiometryValue[1] = "";
					stoichiometryValue[2] = this.metaboliteCompartments.get(productID);
					
					//verify if compartments of metabolites in reaction matches
//					String compID = productID.split("_")[2];
					String compID = null;
					for(String component : productID.split("_"))
						if(this.compartments.keySet().contains(component))
							compID=component;
					
					if(compID!=null && !compID.isEmpty()){
						
						String compartmentName = this.compartments.get(compID);
						
//						if(compartmentName.contains("_"))
//							compartmentName = compartmentName.split("_")[0];
						
						if(!compartmentName.equals(this.metaboliteCompartments.get(productID))){
							if(reactionContainer.getNotes()==null || reactionContainer.getNotes().isEmpty())
								reactionContainer.setNotes("verify products compartments for this reaction");
						}
					}
					
					productsStoichiometry.put(productID, stoichiometryValue);
				}
			}
			
			reactionContainer.setProductsStoichiometry(productsStoichiometry);
			
			//PATHWAYS
			if(drains.contains(idReaction)){
				
				Set<String> pathway = new HashSet<>();
				pathway.add("D0001");
				reactionContainer.setPathways(pathway);
				
				Map<String, String> pathwaysMap = new HashMap<>();
				pathwaysMap.put("D0001", "Drains pathway");
				
				reactionContainer.setPathwaysMap(pathwaysMap);
			}
			
			else if(transportReactions.contains(idReaction)){
				
				Set<String> pathway = new HashSet<>();
				pathway.add("T0001");
				reactionContainer.setPathways(pathway);
				
				Map<String, String> pathwaysMap = new HashMap<>();
				pathwaysMap.put("T0001", "Transporters pathway");
				
				reactionContainer.setPathwaysMap(pathwaysMap);
			}
				
			else if(keggPathwaysData.existsReactionIDinKeggPathway(externalID)){
				
				Set<String> pathways = new HashSet<String>(keggPathwaysData.getReactionPathways(externalID));
				reactionContainer.setPathways(pathways);
				
				Map<String, String> pathwaysMap = new HashMap<>();
				
				for(String pathway : pathways)
					pathwaysMap.put(pathway, keggPathwaysData.getPathwayName(pathway));
				
				reactionContainer.setPathwaysMap(pathwaysMap);
			}

			this.resultReactions.add(reactionContainer);
		}
	}
	

	/**
	 * 
	 */
	public void readEnzymes(){

		//		ENZYME CONTAINER
		//		private String entry;	X
		//		private String name;	X
		//		private List<String> names;	X
		//		private	List<String> dblinks;
		//		private String enzyme_class;	X
		//		private List<String> orthologues;
		//		private	List<String> cofactors;
		//		private List<String> reactions;	X
		//		private	Map<String, String> pathways;
		//		private	List<String> genes;

		Map<String, Set<String>> ecNumbers = cont.getECNumbers();
		
		for(String ecNumber : ecNumbers.keySet()){
			
			EnzymeContainer enzymeContainer = new EnzymeContainer(ecNumber);
			
			Set<String> reactionsSet = ecNumbers.get(ecNumber);
			List<String> reactions = new ArrayList<>();
			List<String> names = new ArrayList<>();
			
			for(String reaction : reactionsSet){
				
				String reactionID = "";
				
				if(reaction.contains("_")){
					
//					String[] reacNameSplit = reaction.toLowerCase().split("_");
					String[] reacNameSplit = reaction.split("_");
					
					if(reaction.startsWith("R"))
						reactionID = reacNameSplit[1];
					else
						reactionID = reacNameSplit[0].concat("_").concat(reacNameSplit[1]);
				}
				else
					reactionID = reaction;
				
				reactions.add(reactionID);
				
				String name = "";
				
				if(this.reactionsData!=null && this.reactionsData.existsReactionID(reactionID))
					name = this.reactionsData.getEnzymeName(reactionID);
				else if(this.enzymeNames.containsKey(reactionID))
					name = this.enzymeNames.get(reactionID);
				
				if(name.length()>10)
					names.add(name);
			}
			
			enzymeContainer.setReactions(reactions);

			if(names.size()==1){
				enzymeContainer.setName(names.get(0));
				enzymeContainer.setNames(names);
			}
			else
				enzymeContainer.setNames(names);
			
			enzymeContainer.setEnzyme_class(ecNumber);

			if(enzymeContainer.getEntryID() != null)
				this.resultEnzymes.add(enzymeContainer);
		}
	}
	
	
	/**
	 * 
	 */
	public void readPathways() {
		
//		PATHWAYS HIERARCHY CONTAINER
//		private Map<String, List<String[]>> pathways_hierarchy;	X
//		private String super_pathway;	X
		
		Map<String , List<String>> superPathways = keggPathwaysData.getSuperPathways();
		
		for(String superPathway : superPathways.keySet()){
			
			PathwaysHierarchyContainer pathwayHierarchyContainer = new PathwaysHierarchyContainer(superPathway);
			
			List<String> intermediaryPathways = superPathways.get(superPathway);
			Map<String, List<String[]>> allPathwaysHierarchy = keggPathwaysData.getPathwaysHierarchy();
			
			Map<String, List<String[]>> pathwaysHierarchy = new HashMap<>();
			
			for(String intermediaryPathway : intermediaryPathways){
				
				List<String[]> pathways = allPathwaysHierarchy.get(intermediaryPathway);
				pathwaysHierarchy.put(intermediaryPathway, pathways);
			}
			
			pathwayHierarchyContainer.setPathways_hierarchy(pathwaysHierarchy);
			
			this.resultPathwaysHierarchy.add(pathwayHierarchyContainer);
		}
	}


	/**
	 * 
	 */
	public void readCompartments(){

		//		COMPARTMENT CONTAINER
		//		private String compartmentID;	X
		//		private String name;	X
		//		private String abbreviation;	X

		Map<String, CompartmentCI> compartments = cont.getCompartments();

		for(String compartmentID : compartments.keySet()){
			
			CompartmentCI compartment = compartments.get(compartmentID);
			String compartmentName = compartment.getName().split("_")[0];

			CompartmentContainer compartmentContainer = new CompartmentContainer(compartmentID, compartmentName, compartmentID);

			this.compartments.put(compartmentID, compartmentName);
			this.resultCompartments.add(compartmentContainer);
		}
	}

	
	
	/**
	 * @param sbmlMetaboliteID
	 * @return
	 */
	public static String processSBMLGeneID(String sbmlGeneID) {
		
		String processedGeneID = sbmlGeneID;
		
		if(processedGeneID.matches("^[Gg]_.+"))
			processedGeneID = processedGeneID.replaceAll("^[Gg]_","");
		
		
		
		return processedGeneID;
	}
	
	/**
	 * @param sbmlMetaboliteID
	 * @return
	 */
	public static String processSBMLMetaboliteID(String sbmlMetaboliteID) {
		
		String processedMetID = sbmlMetaboliteID;
		
		if(processedMetID.matches("^[Mm]_.+"))
			processedMetID = processedMetID.replaceAll("^[Mm]_","");
		
		String[] metIDSplited = sbmlMetaboliteID.split("_");
		
		for(String splitedId : metIDSplited){
			
			if(splitedId.matches("cpd\\d{5}"))
				processedMetID = splitedId;
		}
		
		if(processedMetID.matches(".+[^_]_\\w{1,2}$"))
			processedMetID = processedMetID.substring(0, processedMetID.lastIndexOf("_"));		
		
		return processedMetID;
	}

	
	/**
	 * @return
	 */
	public ConcurrentLinkedQueue<GeneContainer> getResultGenes(){
		return resultGenes;
	}

	/**
	 * @return
	 */
	public ConcurrentLinkedQueue<MetaboliteContainer> getResultMetasbolites(){
		return resultMetabolites;
	}

	/**
	 * @return
	 */
	public ConcurrentLinkedQueue<ReactionContainer> getResultReactions(){
		return resultReactions;
	}

	/**
	 * @return
	 */
	public ConcurrentLinkedQueue<EnzymeContainer> getResultEnzymes(){
		return resultEnzymes;
	}

	/**
	 * @return
	 */
	public ConcurrentLinkedQueue<CompartmentContainer> getResultCompartments(){
		return resultCompartments;
	}
	
	/**
	 * @return
	 */
	public ConcurrentLinkedQueue<PathwaysHierarchyContainer> getResultPathwaysHierarchy(){
		return resultPathwaysHierarchy;
	}
}
