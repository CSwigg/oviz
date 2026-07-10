"""Self-contained USDZ writer used by the Oviz AR snapshot runtime."""


THREEJS_AR_USDZ_RUNTIME_JS = r"""
      // This exporter follows the structure used by the MIT-licensed Three.js
      // USDZExporter, but keeps only the primitives Oviz needs and includes a
      // store-only ZIP writer so AR export has no runtime CDN dependency.
      const OVIZ_USDZ_PRECISION = 7;

      function ovizUsdZNumber(value) {
        const number = Number(value);
        return (Number.isFinite(number) ? number : 0.0).toPrecision(OVIZ_USDZ_PRECISION);
      }

      function ovizUsdZHeader(defaultPrim = "Root") {
        return `#usda 1.0
(
  customLayerData = {
    string creator = "Oviz AR USDZ Exporter"
  }
  defaultPrim = "${defaultPrim}"
  metersPerUnit = 1
  upAxis = "Y"
)

`;
      }

      function ovizUsdZColor(color, fallback = { r: 1, g: 1, b: 1 }) {
        const source = color && typeof color === "object" ? color : fallback;
        return `(${ovizUsdZNumber(source.r)}, ${ovizUsdZNumber(source.g)}, ${ovizUsdZNumber(source.b)})`;
      }

      function ovizUsdZMatrix(matrix) {
        const values = matrix && matrix.elements ? matrix.elements : [
          1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1,
        ];
        const row = (offset) => `(${[
          values[offset], values[offset + 1], values[offset + 2], values[offset + 3],
        ].map(ovizUsdZNumber).join(", ")})`;
        return `(${row(0)}, ${row(4)}, ${row(8)}, ${row(12)})`;
      }

      function ovizUsdZVector3Array(attribute, count) {
        if (!attribute) {
          return Array(Math.max(0, Number(count) || 0)).fill("(0, 0, 0)").join(", ");
        }
        const values = [];
        for (let index = 0; index < attribute.count; index += 1) {
          values.push(`(${ovizUsdZNumber(attribute.getX(index))}, ${ovizUsdZNumber(attribute.getY(index))}, ${ovizUsdZNumber(attribute.getZ(index))})`);
        }
        return values.join(", ");
      }

      function ovizUsdZVector2Array(attribute) {
        const values = [];
        for (let index = 0; index < attribute.count; index += 1) {
          values.push(`(${ovizUsdZNumber(attribute.getX(index))}, ${ovizUsdZNumber(1.0 - attribute.getY(index))})`);
        }
        return values.join(", ");
      }

      function ovizUsdZFaceVertexIndices(geometry) {
        const values = [];
        if (geometry.index) {
          for (let index = 0; index < geometry.index.count; index += 1) {
            values.push(Math.round(geometry.index.getX(index)));
          }
        } else {
          const count = geometry.attributes.position.count;
          for (let index = 0; index < count; index += 1) {
            values.push(index);
          }
        }
        return values.join(", ");
      }

      function ovizUsdZGeometryFile(geometry) {
        if (!geometry || !geometry.attributes || !geometry.attributes.position) {
          throw new Error("AR geometry has no position attribute.");
        }
        if (!geometry.attributes.normal && typeof geometry.computeVertexNormals === "function") {
          geometry.computeVertexNormals();
        }
        const position = geometry.attributes.position;
        const normal = geometry.attributes.normal;
        const uv = geometry.attributes.uv;
        const color = geometry.attributes.color;
        const faceIndexCount = geometry.index ? geometry.index.count : position.count;
        const triangleCount = Math.floor(faceIndexCount / 3);
        const faceCounts = Array(triangleCount).fill("3").join(", ");
        const primvars = [];
        if (uv) {
          primvars.push(`    texCoord2f[] primvars:st = [${ovizUsdZVector2Array(uv)}] (
      interpolation = "vertex"
    )`);
        }
        if (color) {
          primvars.push(`    color3f[] primvars:displayColor = [${ovizUsdZVector3Array(color, color.count)}] (
      interpolation = "vertex"
    )`);
        }
        const doubleSided = geometry.userData && geometry.userData.ovizArDoubleSided
          ? "    bool doubleSided = 1\n"
          : "";
        return `${ovizUsdZHeader("Geometry")}def Xform "Geometry"
{
  def Mesh "Geometry"
  {
    int[] faceVertexCounts = [${faceCounts}]
    int[] faceVertexIndices = [${ovizUsdZFaceVertexIndices(geometry)}]
    normal3f[] normals = [${ovizUsdZVector3Array(normal, position.count)}] (
      interpolation = "vertex"
    )
    point3f[] points = [${ovizUsdZVector3Array(position, position.count)}]
${doubleSided}${primvars.length ? `${primvars.join("\n")}\n` : ""}    uniform token subdivisionScheme = "none"
  }
}
`;
      }

      function ovizUsdZTextureKey(texture) {
        return `${Math.round(Number(texture && texture.id) || 0)}_${texture && texture.flipY === false ? "n" : "f"}`;
      }

      function ovizUsdZTextureFormat(texture) {
        const requested = String(texture && texture.userData && texture.userData.ovizArTextureFormat || "png").toLowerCase();
        return requested === "jpeg" || requested === "jpg" ? "jpg" : "png";
      }

      function ovizUsdZMaterial(material, textures) {
        const materialId = Math.round(Number(material && material.id) || 0);
        const color = material && material.color ? material.color : { r: 1, g: 1, b: 1 };
        const emissive = material && material.emissive ? material.emissive : { r: 0, g: 0, b: 0 };
        const emissiveIntensity = Number.isFinite(Number(material && material.emissiveIntensity))
          ? Number(material.emissiveIntensity)
          : 1.0;
        const emissiveColor = {
          r: Number(emissive.r || 0) * emissiveIntensity,
          g: Number(emissive.g || 0) * emissiveIntensity,
          b: Number(emissive.b || 0) * emissiveIntensity,
        };
        const opacity = Math.min(Math.max(Number(material && material.opacity), 0.0), 1.0);
        const safeOpacity = Number.isFinite(opacity) ? opacity : 1.0;
        const roughness = Number.isFinite(Number(material && material.roughness)) ? Number(material.roughness) : 0.72;
        const metalness = Number.isFinite(Number(material && material.metalness)) ? Number(material.metalness) : 0.0;
        const map = material && material.map ? material.map : null;
        const inputs = [];
        const textureNodes = [];
        if (map) {
          const textureKey = ovizUsdZTextureKey(map);
          const extension = ovizUsdZTextureFormat(map);
          const textureName = `Texture_${textureKey}`;
          textures.set(textureKey, { texture: map, extension });
          inputs.push(`      color3f inputs:diffuseColor.connect = </Materials/Material_${materialId}/${textureName}.outputs:rgb>`);
          if (material.transparent === true || safeOpacity < 1.0) {
            inputs.push(`      float inputs:opacity.connect = </Materials/Material_${materialId}/${textureName}.outputs:a>`);
          }
          textureNodes.push(`
    def Shader "PrimvarReader_diffuse"
    {
      uniform token info:id = "UsdPrimvarReader_float2"
      float2 inputs:fallback = (0.0, 0.0)
      string inputs:varname = "st"
      float2 outputs:result
    }

    def Shader "${textureName}"
    {
      uniform token info:id = "UsdUVTexture"
      asset inputs:file = @textures/${textureName}.${extension}@
      float2 inputs:st.connect = </Materials/Material_${materialId}/PrimvarReader_diffuse.outputs:result>
      token inputs:sourceColorSpace = "sRGB"
      token inputs:wrapS = "clamp"
      token inputs:wrapT = "clamp"
      float outputs:r
      float outputs:g
      float outputs:b
      float outputs:a
      float3 outputs:rgb
    }`);
        } else {
          inputs.push(`      color3f inputs:diffuseColor = ${ovizUsdZColor(color)}`);
          inputs.push(`      float inputs:opacity = ${ovizUsdZNumber(safeOpacity)}`);
        }
        if (Math.max(emissiveColor.r, emissiveColor.g, emissiveColor.b) > 0) {
          inputs.push(`      color3f inputs:emissiveColor = ${ovizUsdZColor(emissiveColor)}`);
        }
        inputs.push(`      float inputs:roughness = ${ovizUsdZNumber(roughness)}`);
        inputs.push(`      float inputs:metallic = ${ovizUsdZNumber(metalness)}`);
        inputs.push("      int inputs:useSpecularWorkflow = 0");
        return `  def Material "Material_${materialId}"
  {
    def Shader "PreviewSurface"
    {
      uniform token info:id = "UsdPreviewSurface"
${inputs.join("\n")}
      token outputs:surface
    }

    token outputs:surface.connect = </Materials/Material_${materialId}/PreviewSurface.outputs:surface>
${textureNodes.join("\n")}
  }
`;
      }

      function ovizUsdZXform(object, geometry, material) {
        const objectId = Math.round(Number(object && object.id) || 0);
        const geometryId = Math.round(Number(geometry && geometry.id) || 0);
        const materialId = Math.round(Number(material && material.id) || 0);
        return `      def Xform "Object_${objectId}" (
        prepend references = @./geometries/Geometry_${geometryId}.usda@</Geometry>
        prepend apiSchemas = ["MaterialBindingAPI"]
      )
      {
        matrix4d xformOp:transform = ${ovizUsdZMatrix(object.matrixWorld)}
        uniform token[] xformOpOrder = ["xformOp:transform"]
        rel material:binding = </Materials/Material_${materialId}>
      }
`;
      }

      function ovizUsdZImageCanvas(image, flipY, maxTextureSize) {
        if (!image || !(Number(image.width) > 0) || !(Number(image.height) > 0)) {
          throw new Error("AR texture has no usable image data.");
        }
        const scale = Math.min(1.0, Math.max(1, Number(maxTextureSize) || 2048) / Math.max(image.width, image.height));
        const canvasEl = document.createElement("canvas");
        canvasEl.width = Math.max(1, Math.round(image.width * scale));
        canvasEl.height = Math.max(1, Math.round(image.height * scale));
        const context = canvasEl.getContext("2d", { alpha: true });
        if (!context) {
          throw new Error("AR texture canvas is unavailable.");
        }
        if (flipY === true) {
          context.translate(0, canvasEl.height);
          context.scale(1, -1);
        }
        context.drawImage(image, 0, 0, canvasEl.width, canvasEl.height);
        return canvasEl;
      }

      function ovizUsdZCanvasBytes(canvasEl, mimeType, quality) {
        return new Promise((resolve, reject) => {
          canvasEl.toBlob(async (blob) => {
            if (!blob) {
              reject(new Error("AR texture encoding failed."));
              return;
            }
            resolve(new Uint8Array(await blob.arrayBuffer()));
          }, mimeType, quality);
        });
      }

      let ovizUsdZCrcTable = null;

      function ovizUsdZCrc32(bytes) {
        if (!ovizUsdZCrcTable) {
          ovizUsdZCrcTable = new Uint32Array(256);
          for (let index = 0; index < 256; index += 1) {
            let value = index;
            for (let bit = 0; bit < 8; bit += 1) {
              value = (value & 1) ? (0xedb88320 ^ (value >>> 1)) : (value >>> 1);
            }
            ovizUsdZCrcTable[index] = value >>> 0;
          }
        }
        let crc = 0xffffffff;
        for (let index = 0; index < bytes.length; index += 1) {
          crc = ovizUsdZCrcTable[(crc ^ bytes[index]) & 0xff] ^ (crc >>> 8);
        }
        return (crc ^ 0xffffffff) >>> 0;
      }

      function ovizUsdZWriteStoredZip(files) {
        if (!Array.isArray(files) || !files.length) {
          throw new Error("USDZ package has no files.");
        }
        const encoder = new TextEncoder();
        const entries = [];
        let localOffset = 0;
        files.forEach((file) => {
          const nameBytes = encoder.encode(String(file.name || ""));
          const dataBytes = file.data instanceof Uint8Array ? file.data : new Uint8Array(file.data || 0);
          const baseDataOffset = localOffset + 30 + nameBytes.length;
          let extraLength = 0;
          if ((baseDataOffset & 63) !== 0) {
            extraLength = (64 - ((baseDataOffset + 4) & 63)) & 63;
            extraLength += 4;
          }
          const dataOffset = baseDataOffset + extraLength;
          if ((dataOffset & 63) !== 0) {
            throw new Error("USDZ ZIP alignment failed.");
          }
          const entry = {
            nameBytes,
            dataBytes,
            extraLength,
            localOffset,
            crc: ovizUsdZCrc32(dataBytes),
          };
          entries.push(entry);
          localOffset = dataOffset + dataBytes.length;
        });
        const centralOffset = localOffset;
        const centralSize = entries.reduce((sum, entry) => sum + 46 + entry.nameBytes.length, 0);
        const output = new Uint8Array(centralOffset + centralSize + 22);
        const view = new DataView(output.buffer);
        const set16 = (offset, value) => view.setUint16(offset, value, true);
        const set32 = (offset, value) => view.setUint32(offset, value >>> 0, true);
        entries.forEach((entry) => {
          let offset = entry.localOffset;
          set32(offset, 0x04034b50); offset += 4;
          set16(offset, 20); offset += 2;
          set16(offset, 0x0800); offset += 2;
          set16(offset, 0); offset += 2;
          set16(offset, 0); offset += 2;
          set16(offset, 33); offset += 2;
          set32(offset, entry.crc); offset += 4;
          set32(offset, entry.dataBytes.length); offset += 4;
          set32(offset, entry.dataBytes.length); offset += 4;
          set16(offset, entry.nameBytes.length); offset += 2;
          set16(offset, entry.extraLength); offset += 2;
          output.set(entry.nameBytes, offset); offset += entry.nameBytes.length;
          if (entry.extraLength) {
            set16(offset, 0x4f56);
            set16(offset + 2, entry.extraLength - 4);
            offset += entry.extraLength;
          }
          output.set(entry.dataBytes, offset);
        });
        let centralCursor = centralOffset;
        entries.forEach((entry) => {
          let offset = centralCursor;
          set32(offset, 0x02014b50); offset += 4;
          set16(offset, 20); offset += 2;
          set16(offset, 20); offset += 2;
          set16(offset, 0x0800); offset += 2;
          set16(offset, 0); offset += 2;
          set16(offset, 0); offset += 2;
          set16(offset, 33); offset += 2;
          set32(offset, entry.crc); offset += 4;
          set32(offset, entry.dataBytes.length); offset += 4;
          set32(offset, entry.dataBytes.length); offset += 4;
          set16(offset, entry.nameBytes.length); offset += 2;
          set16(offset, 0); offset += 2;
          set16(offset, 0); offset += 2;
          set16(offset, 0); offset += 2;
          set16(offset, 0); offset += 2;
          set32(offset, 0); offset += 4;
          set32(offset, entry.localOffset); offset += 4;
          output.set(entry.nameBytes, offset);
          centralCursor += 46 + entry.nameBytes.length;
        });
        let end = centralOffset + centralSize;
        set32(end, 0x06054b50); end += 4;
        set16(end, 0); end += 2;
        set16(end, 0); end += 2;
        set16(end, entries.length); end += 2;
        set16(end, entries.length); end += 2;
        set32(end, centralSize); end += 4;
        set32(end, centralOffset); end += 4;
        set16(end, 0);
        return output;
      }

      class OvizUSDZExporter {
        async parse(sceneAr, options = {}) {
          const config = Object.assign({ maxTextureSize: 4096, jpegQuality: 0.88 }, options || {});
          if (!sceneAr || typeof sceneAr.traverseVisible !== "function") {
            throw new Error("AR scene is unavailable.");
          }
          if (typeof sceneAr.updateMatrixWorld === "function") {
            sceneAr.updateMatrixWorld(true);
          }
          const geometryById = new Map();
          const materialById = new Map();
          const xforms = [];
          sceneAr.traverseVisible((object) => {
            if (!object || !object.isMesh || !object.geometry || !object.material) {
              return;
            }
            const material = Array.isArray(object.material) ? object.material[0] : object.material;
            if (!material) {
              return;
            }
            geometryById.set(Math.round(Number(object.geometry.id) || 0), object.geometry);
            materialById.set(Math.round(Number(material.id) || 0), material);
            xforms.push(ovizUsdZXform(object, object.geometry, material));
          });
          if (!xforms.length) {
            throw new Error("AR scene contains no exportable meshes.");
          }
          const textures = new Map();
          const materialsText = Array.from(materialById.values())
            .map((material) => ovizUsdZMaterial(material, textures))
            .join("\n");
          const sceneText = `${ovizUsdZHeader()}def Xform "Root"
{
  def Scope "Scenes" (
    kind = "sceneLibrary"
  )
  {
    def Xform "Scene" (
      customData = {
        bool preliminary_collidesWithEnvironment = 0
        string sceneName = "Oviz AR Snapshot"
      }
      sceneName = "Oviz AR Snapshot"
    )
    {
      token preliminary:anchoring:type = "plane"
      token preliminary:planeAnchoring:alignment = "horizontal"
${xforms.join("\n")}    }
  }
}

def Scope "Materials"
{
${materialsText}}
`;
          const encoder = new TextEncoder();
          const files = [{ name: "model.usda", data: encoder.encode(sceneText) }];
          Array.from(geometryById.entries()).forEach(([geometryId, geometry]) => {
            files.push({
              name: `geometries/Geometry_${geometryId}.usda`,
              data: encoder.encode(ovizUsdZGeometryFile(geometry)),
            });
          });
          for (const [textureKey, textureEntry] of textures.entries()) {
            const canvasEl = ovizUsdZImageCanvas(
              textureEntry.texture.image,
              textureEntry.texture.flipY !== false,
              config.maxTextureSize
            );
            const isJpeg = textureEntry.extension === "jpg";
            files.push({
              name: `textures/Texture_${textureKey}.${textureEntry.extension}`,
              data: await ovizUsdZCanvasBytes(
                canvasEl,
                isJpeg ? "image/jpeg" : "image/png",
                isJpeg ? config.jpegQuality : 1.0
              ),
            });
          }
          const bytes = ovizUsdZWriteStoredZip(files);
          return bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
        }
      }
""".strip()
