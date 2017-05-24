#coding: utf-8
import requests
import xml.etree.ElementTree as ET
from io import BytesIO
import PIL.Image
from PIL import Image, ImageDraw,ImageFilter
import imageio
import numpy as np

def get_xml_element_of_kegg(content_xml):
    organism_id = content_xml.get('org').upper()
    #pathway_id = content_xml.get('number')
    pathway_title = content_xml.get('title')
    print("get pathway: %s\n"%pathway_title)
    for entry_node in content_xml.findall('entry'):
        entity_ids = tuple(entry_node.get("name").upper().split())
        entity_type = entry_node.get("type")
        
        graphic_element_node = entry_node.find("graphics")
        graphic_element_type = graphic_element_node.get("type")
        
        if graphic_element_type in ("rectangle", "circle"):
            w = int(graphic_element_node.get("width"))
            h = int(graphic_element_node.get("height"))
            x0 = int(graphic_element_node.get("x")) - w/2
            y0 = int(graphic_element_node.get("y")) - h/2
            x1, y1 = x0 + w, y0 + h
            graphic_element = (graphic_element_type, x0, y0, x1, y1)
            
        else:
            print("ignoring grahic element '%s'\n" % graphic_element_type)
            continue
        if entity_type == "gene":
            for entity_id in entity_ids:
                if (not entity_id.startswith(organism_id)):
                    print("invalid gene identifier: %s\n" % entity_id)
                yield(entity_id[len(organism_id) + 1:], entity_type, graphic_element)
        elif entity_type == 'compound':
            for entity_id in entity_ids:
                namespace = None
                for ns in ("CPD", "GL"):
                    if entity_id.startswith(ns + ":"):
                        namespace = ns
                        break
                if namespace == None:
                    print("unkown entity type: %s\n" % entity_id)
                    continue
                yield (entity_id[len(ns) + 1:], entity_type, graphic_element)
        else:
            print("ignoring entity '%s'\n" % entity_type)
            continue

def get_req_xml_kegg(path_id):
    url = "http://rest.kegg.jp/get/%s/kgml" % (path_id)
    
    content_req = requests.get(url)
    content_xml = ET.parse(BytesIO(content_req.content)).getroot()
    
    image_url = content_xml.get("image")
    image_req = requests.get(image_url)
    pathway_image = PIL.Image.open(BytesIO(image_req.content))
    return(pathway_image, content_xml)

def parse_entitys_dict(entities):
    entity_graphic_dict = {}
    graphic_entity_dict = {}
    for (entity_id, entity_type, graphic_element) in entities:
        entity_graphic_dict\
            .setdefault((entity_id, entity_type),[])\
            .append(graphic_element)
        #graphic_entity_dict\
            #.setdefault(graphic_element,[])\
            #.append((entity_id, entity_type))
    return(entity_graphic_dict)#, graphic_entity_dict)

def extract_graphic_entity(entry_items, color_list, entity_graphic_dict, etype='gene', label_list=None):
    for i in range(len(entry_items)):
        if (entry_items[i], etype) in entity_graphic_dict.keys():
            for shape_item, x0, x1, y0, y1 in entity_graphic_dict[(entry_items[i], etype)]:
                if label_list:
                    yield (shape_item, x0, x1, y0, y1, color_list[i], entry_items[i], label_list[i])
                else:
                    yield (shape_item, x0, x1, y0, y1, color_list[i], entry_items[i], None)
        else:
            print("entry ignored %s\n" % (entry_items[i]))

def gif_make(movie, pathway_image, graphic_element_list, duration=11):
    
    for j in range(duration):
        print(str(j+1) + ' picture generated!!!')
        pathway_image_frame = pathway_image.copy()
        for i in range(len(graphic_element_list)):
            padding = 1
            (graphic_element_type, x0, y0, x1, y1, colors, geneid, label) = graphic_element_list[i]
            colors = tuple([int(colors[0] + (255-colors[0])/duration*(j+1)), 
                   int(colors[1] + (255-colors[1])/duration*(j+1)),
                   int(colors[2] + (255-colors[2])/duration*(j+1))])
            
            graphic_element_image = Image.new("RGBA",
                (int((x1 - x0)*1 + padding * 2) + 1,
                int((y1 - y0)*1 + padding * 2) + 1),
                color = colors)
#            graphic_element_image_canvas = PIL.ImageDraw(graphic_element_image)
            
            graphic_element_mask = Image.new("L",
                graphic_element_image.size, color = 255)

            graphic_element_mask_canvas = ImageDraw.Draw(
                graphic_element_mask)
            graphic_element_bbox = (
                padding, padding,
                (x1 - x0) + padding,
                (y1 - y0) + padding)
            
            if graphic_element_type == 'circle':
                graphic_element_mask_canvas.ellipse(
                    graphic_element_bbox,
                    outline = 0, fill = 0)
            else:
                graphic_element_mask_canvas.rectangle(
                    graphic_element_bbox,
                    outline = 0, fill = 0)
            graphic_element_mask = graphic_element_mask.filter(\
                    ImageFilter.GaussianBlur(\
                    radius = 15))
            
            graphic_element_mask = graphic_element_mask.point(
                lambda v: (v) * (1 - 0/ 255))
            graphic_element_image.putalpha(graphic_element_mask)
            pathway_image_frame.paste(
                graphic_element_image,
                (int(x0 - padding  - 1),
                 int(y0 - padding  - 1)),
                graphic_element_image)
        if label:
            pathway_image_frame_canvas = ImageDraw.Draw(pathway_image_frame)
            import PIL.ImageFont
            label_font = PIL.ImageFont.truetype(
                font=r'C:\Windows\Fonts\simhei.ttf', size=20)
            
            for p in range(len(graphic_element_list)):
                (graphic_element_type, x0, y0, x1, y1, colors, geneid, label) = graphic_element_list[p]
                label_width, label_height = label_font.getsize(label)
                label_x = x0 + (x1 - x0)/2 - label_width/2
                label_y = y0 + (y1 - y0)/2 - label_height/2
                pathway_image_frame_canvas.text(
                    (label_x, label_y),
                    label,font=label_font,
                    fill=colors)
        movie.append_data(np.array(pathway_image_frame))
        
def make_gif_kegg(path_id, entry_items, color_list, label_list=None, etype ='gene'):
    
    pathway_image, content_xml = get_req_xml_kegg(path_id)
    entity_graphic_dict = parse_entitys_dict(get_xml_element_of_kegg(content_xml))

    target_list = [i for i in \
               extract_graphic_entity(entry_items, color_list, entity_graphic_dict, etype, label_list)]

    del entity_graphic_dict
    movie = imageio.save(path_id + '.gif', fps=25)
    gif_make(movie, pathway_image, target_list, duration=13)
    del movie
    

make_gif_kegg(path_id = 'hsa03018', entry_items = ['51013','51690','51691'],\
     color_list = [(255,0,0), (0,0,255),(0,0,255)], \
     label_list=['51013','51690','51691'])
