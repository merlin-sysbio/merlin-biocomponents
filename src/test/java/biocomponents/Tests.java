package biocomponents;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.stream.XMLStreamException;

import org.xml.sax.SAXException;

import com.fasterxml.jackson.databind.ObjectMapper;

import pt.uminho.ceb.biosystems.merlin.biocomponents.io.readers.ContainerBuilder;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.DatabaseAccess;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.MySQLDatabaseAccess;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.ContainerUtils;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.ErrorsException;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.JSBMLReader;
import pt.uminho.ceb.biosystems.mew.biocomponents.validation.io.JSBMLValidationException;

/**
 * @author Oscar
 *
 */
public class Tests {

//	@Test
	public void testMerlinToContainer() throws IOException, Exception {
		
		String database = "r6", organism = "spn_r6";

		Container conn = null;//new Container(new JSBMLReader("C:/Users/Oscar Dias/Desktop/spe.xml",organism));
		
//		System.out.println(conn.getBiomassId());
		
		DatabaseAccess connector = new MySQLDatabaseAccess("odias", "password", "192.168.1.85", 3306, database);
		
		conn = new Container(new ContainerBuilder(connector, database, true, organism, ""));
		
		System.out.println(conn.getBiomassId());
		System.out.println(conn.identifyDrains());
		System.out.println(conn.getMetabolitesExtraInfo());
		
		//for(String drain: conn.getDefaultEC().keySet())
			//System.out.println(drain+" "+conn.getReaction(drain).getName()+" "+conn.getDefaultEC().get(drain));
		
	}
	
	public void testMerlinReader() throws IOException, Exception {

		DatabaseAccess m = new MySQLDatabaseAccess("odias", "password", "192.168.1.85", 3306, "kla_model");
		String modelName = "test";
		boolean compartimentalized = true;
		
		Container container = new Container(new ContainerBuilder(m,modelName,compartimentalized,"kla",""));
//		System.out.println(container.getMetabolitesWithoutPathway().size() + "\t" + container.getMetabolitesWithoutPathway());
//		System.out.println(ContainerUtils.getTransportedMetabolites(container, container.getReactions().keySet()));
//		container.getExternalCompartment().getMetabolitesInCompartmentID();
		container.constructDrains(container.getExternalCompartment().getMetabolitesInCompartmentID(), container.getExternalCompartmentId(), 0.0,10000.0);
		
		Set<String> reactions = ContainerUtils.identyfyReactionWithDeadEnds(container);
		for(String reaction : reactions)
			System.out.println(container.getReaction(reaction).getName()/*+"\t"+ContainerUtils.getReactionToString(container, reaction)*/);
		
		System.out.println(container.getMetabolite("M_00807").getReactionsId());
		System.out.println(container.getDrains());
	}
	
	//	@Test
	public void readSBML() throws FileNotFoundException, XMLStreamException, ErrorsException, IOException, ParserConfigurationException, SAXException, JSBMLValidationException {
		
		JSBMLReader jsbml = new JSBMLReader("D:/Dropbox/PAPERS/merlin-manuscript/Revision_1/pylori_analysis/iIT341_min_medI.xml", "hpy", false);
		System.out.println(jsbml.getGenes());
		
		JSBMLReader jsbml2 = new JSBMLReader("D:/Dropbox/PAPERS/merlin-manuscript/Revision_1/pylori_analysis/iIT341_min_medI.xml", "hpy", false);
		System.out.println(jsbml2.getGenes());
		
		JSBMLReader jsbml3 = new JSBMLReader("D:/Dropbox/PAPERS/merlin-manuscript/Revision_1/pylori_analysis/iIT341_min_medI.xml", "hpy", false);
		System.out.println(jsbml3.getGenes());
	}
	
	
//	@Test
	public void testMiriamLib() {
		
		try {
		     URL url = new URL("http://identifiers.org/rest/collections");
		     ObjectMapper om = new ObjectMapper();
		     List<Object> o = om.readValue(url, List.class);
		     List<String> prefixos = new ArrayList<>();
		     for (Object oo : o) {
		       Map<String, Object> map = (Map) oo;
		       prefixos.add(map.get("prefix").toString());
		     }
		     Collections.sort(prefixos);
		     
		     for(String prefix : prefixos)
		    	 System.out.println(prefix);
		     
//		     System.out.println(o.size());
		   } catch (Exception e) {
		     
		   }
		
	}
	
}
