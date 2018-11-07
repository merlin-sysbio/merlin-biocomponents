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
	private List<String> biomassPathway;
	private ModelSources modelSource;
	private Set<String> metabolitesSet, reactionsSet;
//	private Map<String,Integer> genesIds;
	private Map<String,String> compartments;
	private Map<String,Set<String>> reactionEnzymes;
//	private boolean addSpontaneousReactions;
	private String biomassReaction;


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
		this.reactionsSet = new HashSet<>();
		this.compartments = new HashMap<>();
		
//		if(statement!=null)
//			this.genesIds = ModelAPI.getGeneIds(statement);
		
		this.processBiomass();
		
		this.resultPathwaysHierarchy = new ConcurrentLinkedQueue<>();
		this.transportReactions = new ArrayList<>(cont.getReactionsByType(ReactionTypeEnum.Transport));
		this.drains = new ArrayList<>(cont.getReactionsByType(ReactionTypeEnum.Drain));
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
	private void processBiomass(){
		
		this.biomassReaction = cont.getBiomassId();
		
		if(biomassReaction!=null && !biomassReaction.isEmpty())
			if(!cont.getReaction(biomassReaction).getType().equals(ReactionTypeEnum.Biomass)
					&& biomassReaction.toLowerCase().contains("biomass"))
				cont.getReaction(biomassReaction).setType(ReactionTypeEnum.Biomass);
		
		this.biomassPathway = new ArrayList<>(cont.getReactionsByType(ReactionTypeEnum.Biomass));
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
		
		boolean saveDrains = true;
		if(this.drains!=null && !this.drains.isEmpty())
			saveDrains = false;
		
		
//		Set<String> coreMetabolites = new HashSet<>();
//		if(this.metabolitesData!=null)
//			coreMetabolites = this.metabolitesData.getCoreCompounds();
		
		for(String metID : metabolites.keySet()){
			
			if(!metID.endsWith("_b")){

				MetaboliteCI metabolite = metabolites.get(metID);
				
				String metaboliteID = processSBMLMetaboliteID(metID);//.split("_")[1];
				
				String metaboliteCompartmentID = cont.getMetaboliteCompartments(metID).toArray()[0].toString();
				String compartmentName = this.compartments.get(metaboliteCompartmentID);
//				if(metaboliteCompartmentID!=null && !metaboliteCompartmentID.isEmpty())
//					compartmentName = this.compartments.get(metaboliteCompartmentID);
				
				this.metaboliteCompartments.put(metID, compartmentName);


				if(!this.metabolitesSet.contains(metaboliteID)){
					
					MetaboliteContainer metContainer = new MetaboliteContainer(metaboliteID);
					
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
			else if(saveDrains){
				Set<String> reactionsList = cont.getMetabolite(metID).getReactionsId();
				
				for(String drainID : reactionsList)
					if(!this.drains.contains(drainID))
						this.drains.add(drainID);
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
			
			String reactionID = processSBMLReactionID(idReaction);
			
			ReactionCI reaction = reactions.get(idReaction);
			
			Map <String, StoichiometryValueCI> reactants = reaction.getReactants();
			Map <String, StoichiometryValueCI> products = reaction.getProducts();
			
			if((reactants==null || reactants.isEmpty()) && !this.drains.contains(idReaction))
				this.drains.add(idReaction);
			else if((products==null || products.isEmpty()) && !this.drains.contains(idReaction))
				this.drains.add(idReaction);
			
			String compSuffix = "";
			
			if(reactants!=null && !reactants.isEmpty() && reactants.get(reactants.keySet().iterator().next()).getCompartmentId()!=null 
					&& !reactants.get(reactants.keySet().iterator().next()).getCompartmentId().isEmpty()){
				
				compSuffix = compSuffix.concat("_").concat(reactants.get(reactants.keySet().iterator().next()).getCompartmentId());
			}
			
			else if(products!=null && !products.isEmpty() && products.get(products.keySet().iterator().next()).getCompartmentId()!=null 
					&& !products.get(products.keySet().iterator().next()).getCompartmentId().isEmpty()){
				
				compSuffix = compSuffix.concat("_").concat(products.get(products.keySet().iterator().next()).getCompartmentId());
			}
			
			else if(!reaction.identifyCompartments().isEmpty() && reaction.identifyCompartments()!=null){
				
				compSuffix = compSuffix.concat("_").concat(reaction.identifyCompartments().toArray()[0].toString());
			}

			String reactionIDcompartment = reactionID.concat(compSuffix);
			
			int i=1;
			while(this.reactionsSet.contains(reactionIDcompartment)){
				reactionIDcompartment = reactionIDcompartment.concat("_copy")+i;
				i++;
			}
			
			ReactionContainer reactionContainer = new ReactionContainer(reactionIDcompartment);
			reactionContainer.setReversible(reaction.isReversible());
			
			//REACTION COMPARTMENT
			if(!compSuffix.isEmpty())
				reactionContainer.setLocalisation(cont.getCompartment(compSuffix.substring(1)).getName().split("_")[0]);
			
			String externalID = "";
			
			if(this.reactionsData!=null && this.reactionsData.existsReactionID(reactionID)){
				
				if(this.reactionsData.getReactionAbbreviation(reactionID).contains("_"))
					externalID = this.reactionsData.getReactionAbbreviation(reactionID).split("_")[1];
				else
					externalID = this.reactionsData.getReactionAbbreviation(reactionID);
				
				if(externalID.matches("R\\d{5}")){
					
//					try {
//						Integer.parseInt(externalID.substring(1));	//verify we have an ID number with KeggID format ("R"+ numbers)
					reactionIDcompartment = externalID.concat(compSuffix);
					
					i=1;
					while(this.reactionsSet.contains(reactionIDcompartment)){
						reactionIDcompartment = reactionIDcompartment.concat("_copy")+i;
						i++;
					}
					
					reactionContainer.setEntryID(reactionIDcompartment);
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
			
			this.reactionsSet.add(reactionIDcompartment);
			
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
					Double stoichCoef = reactants.get(reactantID).getStoichiometryValue();
					if(stoichCoef>0)
						stoichCoef = -stoichCoef;
					stoichiometryValue[0] = Double.toString(stoichCoef);
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
					Double stoichCoef = products.get(productID).getStoichiometryValue();
					if(stoichCoef<0)
						stoichCoef = -stoichCoef;
					stoichiometryValue[0] = Double.toString(stoichCoef);
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
			Set<String> pathways = new HashSet<>();
			Map<String, String> pathwaysMap = new HashMap<>();
			
			if(idReaction.equals(this.biomassReaction) || this.biomassPathway.contains(idReaction)){
				
				pathways.add("B0001");
				reactionContainer.setPathways(pathways);
				
				pathwaysMap.put("B0001", "Biomass Pathway");
				
				reactionContainer.setPathwaysMap(pathwaysMap);
			}
			
			else if(this.drains.contains(idReaction)){
				
				pathways.add("D0001");
				reactionContainer.setPathways(pathways);
				
				pathwaysMap.put("D0001", "Drains pathway");
				
				reactionContainer.setPathwaysMap(pathwaysMap);
			}
			
			else if(this.transportReactions.contains(idReaction)){
				
				pathways.add("T0001");
				reactionContainer.setPathways(pathways);
				
				pathwaysMap.put("T0001", "Transporters pathway");
				
				reactionContainer.setPathwaysMap(pathwaysMap);
			}
				
			else if(keggPathwaysData.existsReactionIDinKeggPathway(externalID)){
				
				pathways.addAll(keggPathwaysData.getReactionPathways(externalID));
				reactionContainer.setPathways(pathways);
				
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
		
		while(processedMetID.matches("^[Mm]_.+"))
			processedMetID = processedMetID.replaceAll("^[Mm]_","");
		
		if(processedMetID.contains("_")){

			String[] metIDSplited = processedMetID.split("_");

			for(String splitedId : metIDSplited){

				if(splitedId.matches("cpd\\d{5}"))
					processedMetID = splitedId;
			}

			while(processedMetID.matches(".+[^_]_\\w{1}\\d*$"))
				processedMetID = processedMetID.substring(0, processedMetID.lastIndexOf("_"));	

		}
		
		return processedMetID;
	}
	
	
	
		

	
	/**
	 * @param sbmlReactionID
	 * @return
	 */
	public static String processSBMLReactionID(String sbmlReactionID) {
		
		String processedReactionID = sbmlReactionID;
		
		while(processedReactionID.matches("^[Rr]_.+"))
			processedReactionID = processedReactionID.replaceAll("^[Rr]_","");
		
		if(processedReactionID.contains("_")){

			String[] reactionIDSplited = processedReactionID.split("_");

			for(String splitedId : reactionIDSplited){

					if(splitedId.matches("rxn\\d{5}"))
						processedReactionID = splitedId;
				
					else if(splitedId.matches("R\\d{5}"))
						processedReactionID = splitedId;
					
			}

			while(processedReactionID.matches(".+[^_]_\\w{1}\\d*$"))
				processedReactionID = processedReactionID.substring(0, processedReactionID.lastIndexOf("_"));	

		}
		
		return processedReactionID;
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


	/**
	 * @return the biomassReaction
	 */
	public String getBiomassReaction() {
		return biomassReaction;
	}

}
