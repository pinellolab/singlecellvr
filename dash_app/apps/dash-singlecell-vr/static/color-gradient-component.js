/*global THREE, AFRAME, Utils*/

AFRAME.registerComponent('color-gradient', {
    schema: {
        colors: {type: 'array'},
        maxLabel: {type: 'number', default: 1},
        minLabel: {type: 'number', default: 0},
        medianLabel: {type: 'number', default: 0.5},
        height: {type: 'number'},
        width: {type: 'number'},
        verticalOffset: {type: 'number', default: 0},
    },

    init: function () {
        this.height = this.data.height;
        this.width = this.data.width;
        this.colors = this.data.colors;
        this.verticalOffset = -this.data.verticalOffset;
        const blocks = [];
        const blockHeight = 5 / this.colors.length;
        for (let i = 0; i < this.colors.length; i++) {
          const block = Utils.htmlToElement(`<a-entity id=colorblock-${i} position="0 ${this.verticalOffset + (blockHeight * i)} 0" geometry="primitive: plane; height: ${blockHeight}; width: ${this.width};" material="shader: flat; side:front; transparent: true; alphaTest: 0.5; color: ${this.colors[i]}"></a-entity>`);
          blocks.push(block);
        }
        this.el.setAttribute('geometry', `primitive: plane; height: ${this.height}; width: ${this.width};`);
        this.el.setAttribute('material', `shader: flat; side:front; transparent: true; opacity: 0; alphaTest: 0.5; color: white;`);  
        this.el.append(...blocks);
        const scaleBarX = this.width / 2 + .1;
        const scaleBarYBottom = this.verticalOffset - (blockHeight / 2);
        const scaleBarYTop = this.verticalOffset + (blockHeight * (this.colors.length - 1)) + (blockHeight / 2);
        const scaleBarZ = 0;
        const scaleBar = `<a-entity meshline="lineWidth: 5; path: ${scaleBarX} ${scaleBarYBottom} ${scaleBarZ}, ${scaleBarX} ${scaleBarYTop} ${scaleBarZ}; color: black"></a-entity>`;
        this.el.append(Utils.htmlToElement(scaleBar));
        const maxLabelText = `<a-text width="7" value="${this.data.maxLabel}" position="${scaleBarX + .1} ${scaleBarYTop} ${scaleBarZ}" color="black"></a-text>`;
        const minLabelText = `<a-text width="7" value="${this.data.minLabel}" position="${scaleBarX + .1} ${scaleBarYBottom} ${scaleBarZ}" color="black"></a-text>`;
        const medianLabelText = `<a-text width="7" value="${this.data.medianLabel}" position="${scaleBarX + .1} ${(scaleBarYTop + scaleBarYBottom) / 2} ${scaleBarZ}" color="black"></a-text>`;
        this.el.append(Utils.htmlToElement(maxLabelText), Utils.htmlToElement(minLabelText), Utils.htmlToElement(medianLabelText));
    }
});